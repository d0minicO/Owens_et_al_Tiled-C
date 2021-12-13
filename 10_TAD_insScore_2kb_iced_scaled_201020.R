### TILED ANALYSIS ###

library(tidyverse)
library(tools)
library(cowplot)
library(magrittr)
library(RColorBrewer)
library(broom)
library(plotrix)
library(patchwork)

library(FSA)
library(ggpubr)

#### taking iced and scaled interaction matrices calculated by Marieke's script 
# from individual samples

# calculating TAD insulation scores (with replicate info as SD)
# for each bin that is within the TAD
# quantify total reporter interactions that are within the TAD
# or outside the TAD



# created 17 11 2020
# last updated 12 12 2021

###############
### OPTIONS ###
###############

options(max.print=100)


#################
### FUNCTIONS ###
#################


## source an uppererrorbar function from https://stackoverflow.com/questions/27585776/error-bars-for-barplot-only-in-one-direction
source("C:/Users/Dominic/Desktop/Work/Paper/Bioinformatics/RNA-seq/xx_upperErrorBar_plotting.R")

# upper and lower IQR functions for use in stat_summary
iqr_low = function(z) { quantile(z,0.25) }
iqr_hi = function(z) { quantile(z,0.75) }


head_tail <- function(x,description){
  # this function takes a data frame (x) and adds a header containing the "description" string
  # header is added to first row of x
  # tail is added of two empty rows
  
  rbind(c(description, rep("",ncol(x)-1)),
        x,
        rep("",ncol(x)),
        rep("",ncol(x))) %>%
    tibble
  
}



## function from https://stackoverflow.com/questions/15720545/use-stat-summary-to-annotate-plot-with-number-of-observations
## to add text of number of data points ontop of a boxplot 
n_fun <- function(x){
  return(data.frame(y = median(x), label = paste0(length(x))))
}

##############
### INPUTS ###
##############

## public directory
publicDir <- "http://sara.molbiol.ox.ac.uk/public/dowens/CTCF-KO/Tiled/CTCF-KO_virtCapC/"

## 
email <- "dominic.owens@imm.ox.ac.uk"

##
genome <- "mm9"

## the main directory
base = "C:/Users/Dominic/Desktop/Work/Paper/Bioinformatics/Tile-C/CTCF-KO/"


# are we using iced or raw matrices?
matrixType = "iced" # "raw"


####  relative to base folder, or could be specified elsewhere

## the bait file
baitFile <- paste0(base, "baitsRightIDs.txt")

## the fragments genome file
fragFile <- paste0(base, "fragData_2kb.txt")

## the data folder
if (matrixType=="iced"){
  dataFolder <- paste0(base, "iced_matrix/")
} else if (matrixType=="raw") {
  dataFolder <- paste0(base, "matrix/")
} else {
  cat("Warning: matrixType must be set to either \"iced\" or \"raw\"")
}


## tadData
# TAD data containing the location of the
tadFile <- paste0(base, "Runx1_TADs.txt")

## promData
# promoter data including location of bin containing CTCF site deletions
# need this for excluding it later
# format is following cols in tab separated fle
# chr start stop  promName  promBin CTCF_dels_bin
promFile <- paste0(base, "promData.txt")


###############
### OUTPUTS ###
###############

outFolder <- paste0(base, "TAD_insScore_", matrixType, "/")
dir.create(outFolder, showWarnings = F)

variablesFolder <- paste0(base, "savedVariables_", matrixType, "/")
dir.create(variablesFolder, showWarnings = F)

plotFolder <- paste0(base, "TAD_insScore_plots", matrixType, "/")
dir.create(plotFolder, showWarnings = F)

statsFolder = paste0(outFolder, "stats/")
dir.create(statsFolder, showWarnings = F)


##########################
### LOAD THE VARIABLES ###
##########################


cat("Gathering the input files \n")

## load the frag file and make into format useful for joining to later
fragData = data.table::fread(fragFile)
colnames(fragData) <- c("prey_chr", "prey_start", "prey_end", "preyID")
fragData <- 
  fragData %>% 
  mutate("baitID"=preyID, "fragID"=preyID)


## get the TAD data
tadData =  data.table::fread(tadFile)
colnames(tadData) <- c("tad_chr", "prey_start", "prey_end", "tad_name")


## get bin numbers of TAD as vector
tadData %<>%
  left_join(fragData, by="prey_start") %>%
  rename(prey_end=prey_end.x) %>%
  left_join(fragData, by="prey_end") %>%
  select(preyID.x, preyID.y, tad_name)





# get the samples
samples <-
  list.files(path=dataFolder) %>%
  data.frame
colnames(samples) <- "sample"


## make the sample names have . instead of -
samples %<>%
  mutate(sample=gsub("-",".",sample))



##########################
### LOAD THE PROMOTERS ###
##########################
promData=read_tsv(promFile,col_names = c("chr", "start", "stop", "promID", "promName", "delID"))

# get the IDs of the bins with actual CTCF site deletions in to exclude from the counting
# as these would not allow for a fair comparison between genotypes (within genotype would be fine though)
delIDs=promData$delID

# get the proms in a more intuitive order
promData %<>%
  arrange(promName)



#########################
### CHOSE THE SAMPLES ###
#########################


# get the matrix location (sampleFIle) and add to the table
# along with
# tissue types, genotypes, and clone names
samples %<>%
  mutate(sampleFiles=paste0(dataFolder,sample,"/iced/2000/",sample,"_2000_iced.matrix"),
         sample2=sample) %>%
  separate(sample2,into=c("tissue","genotype", "clone"),sep="_")



# get the tissue types
tissues <- levels(factor(word(samples$sample, 1, sep="_")))
# geno types
genotypes <- levels(factor(word(samples$sample, 2, sep="_")))
# exchange the - for . to prevent problems with column naming
genotypes <- gsub("-", ".", genotypes)
# clone names (reps effectively)
clones <- levels(factor(word(samples$sample, 3, sep="_")))





#############################
#### PREVIOUS STEP DATA #####
#############################

# load previously processed long form data
# generated in 8_totalPromContacts_2kb_iced_scaled_200422
d <- readRDS(paste0(variablesFolder, "/Tiled_indiReps_normed_", matrixType, "_longFormat.rds"))



# filter out all the excluded delIDs
for (delID in delIDs){
  d %<>% 
    filter(!baitID %in% delID &
             !preyID %in% delID)
}


## have quite a few zero contacts in the df so remove those as they are meaningless
d %<>%
  filter(contacts>0)



################################################
#### ASSIGN WITHIN OR WITHOUT TAD CONTACTS #####
################################################



total_df = data.frame()
## set up for TAD loop
for (tad in tadData$tad_name){
  
  tadData2=
    tadData %>%
    filter(tad_name==!!tad) %>%
    c %>% unlist
  
  
  # add a column to define which interactions are within the TAD or outside the TAD
  d2 = 
    d %>%
    mutate(TAD=if_else(baitID>=as.numeric(tadData2[1]) & preyID<=as.numeric(tadData2[2]),
                       "IN", "OUT" ))
  
  
  #s=samples$sample[1]
  ## go through each sample one at a time
  
  # set up a main samples df to save into
  full_df = data.frame()
  
  for (s in samples$sample){
    
    sample_filt =
      d2 %>%
      filter(sample==!!s)
    
    
    
    #bait = (tadData2[1]:tadData2[2])[1]
    ## go through each bait, and quantify how many interactions are within the TAD vs without the TAD
    
    # set up a baits df
    bait_df = data.frame()
    
    for (bait in tadData2[1]:tadData2[2]){
      
      bait_filt=
        sample_filt %>%
        filter(baitID==!!bait | preyID==!!bait)
      
      
      
      # quant total INSIDE and OUTSIDE TAD interactions
      ints =
        bait_filt %>%
        group_by(TAD) %>%
        summarise(sum(contacts),.groups = 'drop')
      
      
      
      # construct a df to have outside these loops
      temp_bait_df = 
        data.frame(Bait=bait,
                   IN=as.numeric(ints[1,2]),
                   OUT=as.numeric(ints[2,2]))
      
      
      # get the ratio of IN to OUT interactions
      temp_bait_df %<>%
        mutate(Total=IN+OUT,
               Ratio=IN/Total)
      
      
      # rbind to the main bait_df
      bait_df = rbind(temp_bait_df,bait_df)
      
    } # end of bait loop
    
    # add a sample column to the bait_df
    bait_df %<>%
      mutate(sample=!!s)
    
    
    # combine the sample+bait df to the main df
    full_df = rbind(bait_df,full_df)
    
  } # end of sample loop
  
  
  ### split up the sample column to get better grouping possibilities
  ## add a TAD identified for this TAD
  full_df %<>%
    mutate(sample2=sample,
           tad_name=!!tad) %>%
    separate(sample2,into=c("tissue","genotype", "clone", "XP", "ID"),sep="_")
  
  
  ## only keep complete cases
  full_df = full_df[complete.cases(full_df),]
  
  
  
  ## and join with the total df
  total_df=rbind(full_df,total_df)
  
  
} # end of for TAD loop


### save the variable for next time
#saveRDS(object = total_df, file=paste0(variablesFolder, "_insulationScores_indiBaits.rds"))


### load previous step data
total_df <- readRDS(paste0(variablesFolder, "_insulationScores_indiBaits.rds"))


## get mean insulation scores for each sample
insScores = 
  total_df %>%
  group_by(sample, tissue, genotype, clone, XP, ID, tad_name) %>%
  summarise(mean_ins=mean(Ratio),
            median_ins=median(Ratio))


## make a combined column to find the mean and sd of each tissue_genotype group
insScores %<>%
  mutate(group=paste0(tissue,"_",genotype))


## get a summary table to report in paper
insScores2 = 
  insScores %>%
  ungroup() %>%
  group_by(tissue, genotype, tad_name) %>%
  summarise(mean_median_ins=mean(mean_ins),
            median_median_ins=median(median_ins))



# just the WT summary
insScores2 %>%
  filter(genotype=="WT" & tad_name=="main_TAD")


###################################
########### STATISTICS ############
###################################

# set p adjustment method
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")

meth="holm"


## first need to test for normality to decide whether to use parametric or non-parametric tests

## use exploratory plotting
ggdensity(total_df$Ratio, 
          main = "Density plot of per bin TAD insulation ratios",
          xlab = "per bin TAD insulation ratios")


ggqqplot(total_df$Ratio)

## both look highly non-normal

## doing a shapiro-wilks test of normality
# had to sample as max for func is 5k
shapiro.test(sample(total_df$Ratio,5000))
#W = 0.71561, p-value = 1.486e-12

## data is non-normal and so non-parametric Kruskall Wallis tests should be used



### TEST 1 ###

## WILD TYPE SAMPLES ACROSS DIFFERENTIATION MAIN TAD ONLY

## I know a priori that I am only interested in WT samples now,
# and main TAD
# so do stat test on this smaller subset of the data
# adjusted p values will be given in the Dunn multiple comparisons test

## WT only & P1 only
temp=
  total_df %>% 
  filter(genotype=="WT" & tad_name=="main_TAD")

t1 = kruskal.test(Ratio ~ tissue, data = temp)
print(t1)

#Dunn's Kruskal-Wallis post-hoc test
posthoc1<-dunnTest(Ratio ~ tissue, data=temp, method=meth)
ph1 = 
  print(posthoc1) %>%
  mutate(testGroup=paste0("WT", "_", "main_TAD"),
         KW_pval=tidy(t1)$p.value,
         KW_signif=if_else(tidy(t1)$p.value<0.05, "SIG*", "not"))


## save the table
write.table(ph1, 
            file=paste0(statsFolder, "KW_Dunn_mainTAD_WT_allTissues_insScores.txt"),
            quote=F,
            row.names = F,
            col.names = T)





### TEST 2 ###

## WILD TYPE SAMPLES ACROSS DIFFERENTIATION sub TADs ONLY

## I know a priori that I am only interested in WT samples now,
# and sub TADs
# so do stat test on this smaller subset of the data
# adjusted p values will be given in the Dunn multiple comparisons test
# but as I will be testing the sub-TADs separately, I will combine the two KW tests
# and readjust the pvalues

tads = c("P1-P2_TAD", "P2-3'UTR_TAD")

tests_out = data.frame()
for (tad in tads){
  
  temp = 
    total_df %>%
    filter(tad_name==!!tad & genotype=="WT")
  
  t1 = kruskal.test(Ratio ~ tissue, data = temp)
  #print(t1)
  t1_out = tidy(t1)
  
  t1_out %<>%
    mutate(test_group=tad)
  
  tests_out= rbind.data.frame(tests_out,t1_out)
  
} # end for tad loop


# adjust the KW p values
tests_out$p.adj = p.adjust(tests_out$p.value, method=meth)

# make a new significance column
tests_out %<>% 
  mutate(KW_signif=if_else(p.adj<0.05, "SIG*", "not"))


## both were significant#

# save the KW table
write.table(tests_out, 
            file=paste0(statsFolder, "KW_subTADs_WT_allTissues_insScores.txt"),
            quote=F,
            row.names = F,
            col.names = T)


# now do the Dunn's tests
tests_out2 = data.frame()
for (tad in tads){
  
  
  temp = 
    total_df %>%
    filter(tad_name==!!tad & genotype=="WT")
  
  
  #Dunn's Kruskal-Wallis post-hoc test
  posthoc1<-dunnTest(Ratio ~ tissue, data=temp, method=meth)
  ph1 = 
    print(posthoc1) %>%
    mutate(testGroup=paste0("WT", "_", tad))
  
  
  tests_out2= rbind.data.frame(tests_out2,ph1)
  
} # end for tad loop


## now adjust the pvalues FOR THIS MANY COMPARISONS
phs = tests_out2
phs$p.adj.total = p.adjust(phs$P.unadj, method=meth)

## add a significance column
# and tidy the final table
phs %<>%
  mutate(signif=if_else(p.adj.total<0.05, "SIG*", "not")) %>%
  select(Comparison, testGroup, p.adj.total, signif)



# save the Dunn table
write.table(phs, 
            file=paste0(statsFolder, "KW_Dunn_subTADs_WT_allTissues_insScores.txt"),
            quote=F,
            row.names = F,
            col.names = T)







### TEST 3 ###

## ALL GENOTYPES ACROSS DIFFERENTIATION ALL TADs

## I know a priori that I am only interested in CD41 samples now
# so do stat test on this smaller subset of the data
# adjusted p values will be given in the Dunn multiple comparisons test
# but as I will be testing the sub-TADs separately, I will combine the KW tests
# and readjust the pvalues

tads = c("main_TAD", "P1-P2_TAD", "P2-3'UTR_TAD")

tests_out = data.frame()
for (tad in tads){
  
  temp = 
    total_df %>%
    filter(tad_name==!!tad & tissue=="CD41")
  
  t1 = kruskal.test(Ratio ~ genotype, data = temp)
  #print(t1)
  t1_out = tidy(t1)
  
  t1_out %<>%
    mutate(test_group=paste0(tad,"_CD41"))
  
  tests_out= rbind.data.frame(tests_out,t1_out)
  
} # end for tad loop


# adjust the KW p values
tests_out$p.adj = p.adjust(tests_out$p.value, method=meth)

# make a new significance column
tests_out %<>% 
  mutate(KW_signif=if_else(p.adj<0.05, "SIG*", "not"))


## two were significant#

# save the KW table
write.table(tests_out, 
            file=paste0(statsFolder, "KW_allTADs_allGenos_CD41_insScores.txt"),
            quote=F,
            row.names = F,
            col.names = T)


# now do the Dunn's tests
tests_out2 = data.frame()
for (tad in tads){
  
  
  temp = 
    total_df %>%
    filter(tad_name==!!tad & tissue=="CD41")
  
  
  #Dunn's Kruskal-Wallis post-hoc test
  posthoc1<-dunnTest(Ratio ~ genotype, data=temp, method=meth)
  ph1 = 
    print(posthoc1) %>%
    mutate(testGroup=paste0("CD41", "_", tad))
  
  
  tests_out2= rbind.data.frame(tests_out2,ph1)
  
} # end for tad loop


## now adjust the pvalues FOR THIS MANY COMPARISONS
phs = tests_out2
phs$p.adj.total = p.adjust(phs$P.unadj, method=meth)

## add a significance column
# and tidy the final table
phs %<>%
  mutate(signif=if_else(p.adj.total<0.05, "SIG*", "not")) %>%
  select(Comparison, testGroup, p.adj.total, signif)



# save the Dunn table
write.table(phs, 
            file=paste0(statsFolder, "KW_Dunn_allTADs_allGenos_CD41_insScores.txt"),
            quote=F,
            row.names = F,
            col.names = T)





#################################
########### PLOTTING ############
#################################


# reset the levels to allow nicer plotting
total_df$tissue = factor(total_df$tissue, levels=c("Un", "Flk1", "CD41"))
total_df$genotype = factor(total_df$genotype, levels=c("WT", "P1.CTCF.KO", "P2.CTCF.KO"))


### PLOT 1 ###

# WT just Un and Mes MAIN TAD ONLY

temp = 
  total_df %>%
  filter(genotype=="WT" & (tissue=="Un" | tissue=="Flk1") & tad_name=="main_TAD" )

# check the filtering work a I had some issues
levels(factor(temp$tissue))
levels(factor(temp$genotype))
levels(factor(temp$tad_name))

#check what the median values will be
temp %>%
  group_by(tissue) %>%
  summarise(median(Ratio))


## need to identify which points are considered outliers by a standard boxplot and set the limits to that value
# values returned are lower whisker, lower hinge, median, upper hinge, and upper whisker
stat = boxplot.stats(temp$Ratio)$stats

lims = c(stat[1],stat[5]*1.02)



p=
  ggplot(data=temp, aes(x = tissue, y = Ratio, fill=tissue))+
  geom_boxplot(size=.1,outlier.shape = NA)+
  stat_summary(fun.data = n_fun, geom = "text",color="black",size=1.5,vjust=-.5,aes(y=stat[5]))+
  coord_cartesian(ylim = lims)+
  scale_fill_manual(values=c("#7570b3ff", "#d95f02ff", "#1b9e77ff"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        strip.background =element_rect(fill="white",size=.1),
        axis.ticks = element_line(size=.1),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.1))

ggsave(p, 
       filename = paste0(plotFolder, "Un_Mes_WT_mainTAD_median_IQR.pdf"),
       width=1.85,
       height=1)





### PLOT 2 ###

# WT just Mes and CD41 MAIN TAD ONLY

temp = 
  total_df %>%
  filter(genotype=="WT" & (tissue=="Flk1" | tissue=="CD41") & tad_name=="main_TAD")


# check the filtering work a I had some issues
levels(factor(temp$tissue))
levels(factor(temp$genotype))
levels(factor(temp$tad_name))

#check what the median values will be
temp %>%
  group_by(tissue) %>%
  summarise(median(Ratio))



## need to identify which points are considered outliers by a standard boxplot and set the limits to that value
# values returned are lower whisker, lower hinge, median, upper hinge, and upper whisker
## needs to be done for each group in this case as upper whisker was being cut off when performing boxpot.stats on all data
## need some starting points that are off the scale big and small to allow them to be updated by the real data!!
lwr=1e10
upr=0

tissues = levels(factor(temp$tissue))

for(tiss in tissues){

  temp2 =
    temp %>%
    filter(tissue==!!tiss)
  
  
  stat = boxplot.stats(temp2$Ratio)$stats
  
  lwr_tmp = stat[1]
  upr_tmp = stat[5]
  
  ## update these values in the loop if smaller or bigger
  if(lwr_tmp<lwr){lwr=lwr_tmp}
  if(upr_tmp>upr){upr=upr_tmp}
      
}


lims=c(lwr,upr*1.02)





p=
  ggplot(data=temp, aes(x = tissue, y = Ratio, fill=tissue))+
  geom_boxplot(size=.1,outlier.shape = NA)+
  stat_summary(fun.data = n_fun, geom = "text",color="black",size=1.5,vjust=-.5,aes(y=upr))+
  coord_cartesian(ylim = lims)+
  scale_fill_manual(values=c("#d95f02ff", "#1b9e77ff"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        strip.background =element_rect(fill="white",size=.1),
        axis.ticks = element_line(size=.1),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.1))

ggsave(p, 
       filename = paste0(plotFolder, "Mes_CD41_WT_mainTAD_median_IQR.pdf"),
       width=1.85,
       height=1)





### PLOT 3 ###

# WT just Mes and CD41 sub TADs

temp = 
  total_df %>%
  filter(genotype=="WT" & (tissue=="Flk1" | tissue=="CD41") & tad_name!="main_TAD")

# get a variable for the ylimits
lims=c((min(temp$Ratio)-.03),(max(temp$Ratio)+.01))


# check the filtering work a I had some issues
levels(factor(temp$tissue))
levels(factor(temp$genotype))
levels(factor(temp$tad_name))

#check what the median values will be
temp %>%
  group_by(tissue,tad_name) %>%
  summarise(median(Ratio))


## need to identify which points are considered outliers by a standard boxplot and set the limits to that value
# values returned are lower whisker, lower hinge, median, upper hinge, and upper whisker
stat = boxplot.stats(temp$Ratio)$stats

lims = c(stat[1],stat[5]*1.04)



p=
  ggplot(data=temp, aes(x = tissue, y = Ratio, fill=tissue))+
  geom_boxplot(size=.1,outlier.shape = NA)+
  stat_summary(fun.data = n_fun, geom = "text",color="black",size=1.5,vjust=-.3,aes(y=stat[5]))+
  coord_cartesian(ylim = lims)+
  facet_wrap(~tad_name)+
  scale_fill_manual(values=c("#d95f02ff", "#1b9e77ff"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        strip.background =element_rect(fill="white",size=.1),
        axis.ticks = element_line(size=.1),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.1))

ggsave(p, 
       filename = paste0(plotFolder, "Mes_CD41_WT_subTADs_median_IQR.pdf"),
       width=2.2,
       height=1.25)





### PLOT 4 ###

# just CD41 all TADs across genotypes


## plotting in for loop as I want different lims for each one
tads = c("main_TAD", "P1-P2_TAD", "P2-3'UTR_TAD")

for (tad in tads){
  
  temp = 
    total_df %>%
    filter(tissue=="CD41" & tad_name==!!tad)
  
  ## need to identify which points are considered outliers by a standard boxplot and set the limits to that value
  # values returned are lower whisker, lower hinge, median, upper hinge, and upper whisker
  ## needs to be done for each group in this case as upper whisker was being cut off when performing boxpot.stats on all data
  ## need some starting points that are off the scale big and small to allow them to be updated by the real data!!
  lwr=1e10
  upr=0
  
  genos = levels(factor(temp$genotype))
  
  for(g in genos){
    
    temp2 =
      temp %>%
      filter(genotype==!!g)
    
    
    stat = boxplot.stats(temp2$Ratio)$stats
    
    lwr_tmp = stat[1]
    upr_tmp = stat[5]
    
    ## update these values in the loop if smaller or bigger
    if(lwr_tmp<lwr){lwr=lwr_tmp}
    if(upr_tmp>upr){upr=upr_tmp}
    
  }
  
  
  lims=c(lwr,upr*1.02)
  
  
  
  
  ## put genotypes in a logical order
  temp$genotype = factor(temp$genotype, levels=c("WT","P1.CTCF.KO","P2.CTCF.KO"))
  
  
  # check the filtering work a I had some issues
  levels(factor(temp$tissue))
  levels(factor(temp$genotype))
  levels(factor(temp$tad_name))
  

  p=  
    ggplot(data=temp, aes(x = genotype, y = Ratio, fill=genotype))+
    geom_boxplot(size=.1,outlier.shape = NA)+
    stat_summary(fun.data = n_fun, geom = "text",color="black",size=1.5,vjust=-.4,aes(y=upr))+
    coord_cartesian(ylim = lims)+
    facet_grid(~tad_name, scales="free")+
    scale_fill_manual(values=c("#666666ff", "#66a61eff", "#e7298aff"))+
    ylab(NULL)+
    xlab(NULL)+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          strip.background =element_rect(fill="white",size=.1),
          strip.text = element_text(size=7),
          axis.text = element_text(size=6),
          axis.ticks = element_line(size=.1),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=.1))
  
  ggsave(p, 
         filename = paste0(plotFolder, "CD41_allGenos_",tad,"_boxplot_median_IQR.pdf"),
         width=2.3,
         height=1.2)
  
  
}


