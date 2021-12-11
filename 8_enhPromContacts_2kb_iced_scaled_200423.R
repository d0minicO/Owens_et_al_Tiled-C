### TILED ANALYSIS ###

library(tidyverse)
library(tools)
library(cowplot)
library(magrittr)
library(RColorBrewer)


library(FSA)
library(ggpubr)

#### taking iced and scaled interaction matrices calculated by Marieke's script 
# from individual samples

# calculating enhancer promoter contacts in different groups (with replicate info as SD)


# working 16 11 2020


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


## promData
# promoter data including location of bin containing CTCF site deletions
# need this for excluding it later
# format is following cols in tab separated fle
# chr start stop  promName  promBin CTCF_dels_bin
promFile <- paste0(base, "promData.txt")


## enhData
# had bed file of enhancer cordinates and using bedtools 
# intersected with fragData to get overlapping bins
enhFile <- paste0(base, "enhData_bins_v2.txt")


###############
### OUTPUTS ###
###############

outFolder <- paste0(base, "enhProm_contacts_", matrixType, "/")
dir.create(outFolder, showWarnings = F)

variablesFolder <- paste0(base, "savedVariables_", matrixType, "/")
dir.create(variablesFolder, showWarnings = F)

statsFolder = paste0(outFolder, "stats/")
dir.create(statsFolder, showWarnings = F)

##########################
### LOAD THE VARIABLES ###
##########################


cat("Gathering the input files \n")

## load the bait file, get names of baits, and make into format useful for joining to later
## need to match the wrong DpnII bait IDs
baits = data.table::fread(baitFile)
colnames(baits) <- c("bait_chr", "bait_start", "bait_end", "baitID", "baitName")
baitNames <- baits$baitName
baits <- 
  baits %>% 
  mutate("bait_frag"=paste0(bait_chr,":",bait_start,"-",bait_end)) %>%
  dplyr::select(bait_frag, baitID, baitName)


## load the frag file and make into format useful for joining to later
fragData = data.table::fread(fragFile)
colnames(fragData) <- c("prey_chr", "prey_start", "prey_end", "preyID")
fragData <- 
  fragData %>% 
  mutate("baitID"=preyID, "fragID"=preyID)
#  dplyr::select(prey_frag, preyID)


# get the samples
samples <-
  list.files(path=dataFolder) %>%
  data.frame
colnames(samples) <- "sample"

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


##########################
### LOAD THE ENHANCERS ###
##########################

## load the enhancer data
enhData=read_tsv(enhFile,col_names=c("chr", "start", "stop", "enhID", "chr_enh", "start_enh", "stop_enh", "enhName"))

# prepare the enhancer data table to allow for some enhancers spread over multiple bins
enhancers = levels(factor(enhData$enhName))

enhData_new = data.frame()
for (enhancer in enhancers){
  temp <-
    enhData %>%
    filter(enhName==enhancer) %>%
    mutate(min_enhID=min(enhID),
           max_enhID=max(enhID)) %>%
    mutate(enhID_both=paste0(min_enhID,"_",max_enhID))
  
  enhData_new <- rbind.data.frame(enhData_new,temp)
}

# now filter to make the rows unique
enhData_new %<>%
  distinct(enhID_both, .keep_all = T) %>%
  arrange(enhID)


# check how many bins each enhancer covers
enhData_new %<>%
  mutate(binNum=(max_enhID-min_enhID)+1)

# add a better just numerical enhancer name to use for plotting later
Runx1_ATG=92825890

enhData_new %<>%
  mutate(enh_mid=round((start_enh+stop_enh)/2)) %>%
  mutate(enh_numName=factor(round(((Runx1_ATG-enh_mid)/1e3),digits=1)))

# quickly plot and save this info
binPlot <-
  ggplot(data=enhData_new, aes(x=enh_numName,y=binNum))+
  geom_bar(stat="identity")+
  scale_y_continuous(breaks=(c(0,1,2)),name="Number of bins")+
  scale_x_discrete(name="Enhancer")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

ggsave(plot=binPlot,
       filename=paste0(outFolder,"NumberOfBins_inEachEnhancer.pdf"),
       width=6,
       height=2.5)

# output the enhancer data table
write.table(enhData_new,
            file=paste0(outFolder,"EnhancerData_all.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            na = "",
            sep="\t")

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


############################################
#### GATHER ENHANCER PROMOTER CONTACTS #####
############################################
totals = data.frame()

for (i in 1:nrow(promData)){ # for prom loop
  
  promName=promData$promName[i]
  promID=promData$promID[i]
  
  
  # filter the matrix on just this prom
  prom_filt <-
    d %>% 
    filter(baitID==promID | preyID==promID) %>%
    mutate(promName=!!promName)
  
  # filter out all the excluded delIDs
  for (delID in delIDs){
    prom_filt %<>% 
      filter(!baitID %in% delID &
               !preyID %in% delID)
  }
  
  
  # start the for enhancer loop
  
  for (e in 1:nrow(enhData_new)){
    
    # get the info for just this enhancer
    enh_numName=enhData_new$enh_numName[e]
    enhName=enhData_new$enhName[e]
    min_enhID=enhData_new$min_enhID[e]
    max_enhID=enhData_new$max_enhID[e]
    
    
    # filter on just this enhancers interactions with just this promoter
    # this would fail if enhancers spread over three bins!!!
    enh_filt <-
      prom_filt %>% 
      filter(baitID==min_enhID | baitID==max_enhID |
               preyID==min_enhID | preyID==max_enhID) %>%
      mutate(enh_numName=!!enh_numName,
             enhName=!!enhName)
    
    # add to output df
    totals = rbind.data.frame(totals,enh_filt)
    
  } # end of for enhancer loop
  
} # end of for prom loop





# get total E-P contacts in each group
# start by getting it for each sample
# and then split up the sample into group
sum_EP_contacts <-
  totals %>%
  group_by(sample,promName) %>%
  summarise(sum=sum(contacts)) %>%
  separate(sample,into=c("tissue","genotype", "clone", "XP", "ID"),sep="_") %>%
  mutate(group=paste0(tissue,"_",genotype)) %>%
  ungroup()




###################################
########### STATISTICS ############
###################################


# set p adjustment method
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")

meth="holm"


## first need to test for normality to decide whether to use parametric or non-parametric tests

## use exploratory plotting
ggdensity(totals$contacts, 
          main = "Density plot of E-P contacts",
          xlab = "E-P contacts")


ggqqplot(totals$contacts)

## both look highly non-normal

## doing a shapiro-wilks test of normality
shapiro.test(totals$contacts)
#W = 0.48763, p-value < 2.2e-16

## data is non-normal and so non-parametric Kruskall Wallis tests should be used




### TEST 1 ###

## WILD TYPE SAMPLES ACROSS DIFFERENTIATION

## I know beforehadn that I am only interested in WT samples now,
# so it is okay to do stat test on this smaller subset of the data
# however I am doing tests on P1 and P2 separately
# so doing two separate KS tests
# and then adjusting all the p values together in the Dunn multiple comparisons test

## WT only & P1 only
temp=
  totals %>% 
  filter(genotype=="WT" & promName=="Runx1.P1")

t1 = kruskal.test(contacts ~ tissue, data = temp)
print(t1)

#Dunn's Kruskal-Wallis post-hoc test
posthoc1<-dunnTest(contacts ~ tissue, data=temp, method=meth)
ph1 = 
  print(posthoc1) %>%
  mutate(testGroup=paste0("WT", "_", "Runx1.P1"),
         KS_pval=tidy(t1)$p.value,
         KS_signif=if_else(tidy(t1)$p.value<0.05, "SIG*", "not"))


## WT only & P2 only

temp=
  totals %>% 
  filter(genotype=="WT" & promName=="Runx1.P2")

t2 = kruskal.test(contacts ~ tissue, data = temp)
print(t2)

#Dunn's Kruskal-Wallis post-hoc test
posthoc2<-dunnTest(contacts ~ tissue, data=temp, method=meth)
ph2 = 
  print(posthoc2) %>%
  mutate(testGroup=paste0("WT", "_", "Runx1.P2"),
         KS_pval=tidy(t2)$p.value,
         KS_signif=if_else(tidy(t2)$p.value<0.05, "SIG*", "not"))


## now combine the posthoc tests and adjust the pvalues FOR THIS MANY COMPARISONS
phs = rbind(ph1,ph2)
phs$p.adj.total = p.adjust(phs$P.unadj, method=meth)

## add a significance column
# and tidy the final table
phs %<>%
  mutate(signif=if_else(p.adj.total<0.05, "SIG*", "not")) %>%
  select(Comparison, testGroup, KS_pval, KS_signif, p.adj.total, signif)



## save the table
write.table(phs, 
            file=paste0(statsFolder, "KS_Dunns_adjusted_enhProm_contacts_justWT.txt"),
            quote=F,
            row.names = F,
            col.names = T)



### TEST 2 ###

# ALL GENOTYPES, BOTH PROMOTERS, TOTAL E-Ps


## set a aloop to do these as there will be 6 KS tests

tissues = c("Un", "Flk1", "CD41")
proms = c("Runx1.P1", "Runx1.P2")

tests_out = data.frame()
for (prom in proms){
  
  temp = 
    totals %>%
    filter(promName==!!prom)
  
  for (tiss in tissues){
    
    temp2 =
      temp %>%
      filter(tissue==!!tiss)
    
    
    t1 = kruskal.test(contacts ~ genotype, data = temp2)
    #print(t1)
    t1_out = tidy(t1)
    
    t1_out %<>%
      mutate(test_group=paste0(prom, "_", tiss))
    
    
    tests_out= rbind.data.frame(tests_out,t1_out)
    
  } # end for tissue loop
  
} # end for prom loop


# adjust the KS p values
tests_out$p.adj = p.adjust(tests_out$p.value, method=meth)

# make a new significance column
tests_out %<>% 
  mutate(KS_signif=if_else(p.adj<0.05, "SIG*", "not"))


## none of the KS tests were significant, so not doing any post hoc tests
## save the table
write.table(tests_out, 
            file=paste0(statsFolder, "KS_enhProm_contacts_allGenotypes_allTissues.txt"),
            quote=F,
            row.names = F,
            col.names = T)





### TEST 3 ###

# ALL GENOTYPES, BOTH PROMOTERS, INDIVIDUAL Enhancers


## set a loop to do these as there will be MANY KS tests

tissues = c("Un", "Flk1", "CD41")
proms = c("Runx1.P1", "Runx1.P2")

tests_out = data.frame()
for (prom in proms){
  
  temp = 
    totals %>%
    filter(promName==!!prom)
  
  for (tiss in tissues){
    
    temp2 =
      temp %>%
      filter(tissue==!!tiss)
    
    for (enhancer in enhancers){
      
      temp3 =
        temp2 %>%
        filter(enhName==!!enhancer)
      
      t1 = kruskal.test(contacts ~ genotype, data = temp2)
      #print(t1)
      t1_out = tidy(t1)
      
      t1_out %<>%
        mutate(test_group=paste0(prom, tiss, enhancer, sep="_"))
      
      
      tests_out= rbind.data.frame(tests_out,t1_out)
      
      
    } # end for enhancer loop
    
  } # end for tissue loop
  
} # end for prom loop


# adjust the KS p values
tests_out$p.adj = p.adjust(tests_out$p.value, method=meth)

# make a new significance column
tests_out %<>% 
  mutate(KS_signif=if_else(p.adj<0.05, "SIG*", "not"))


## none of the KS tests were significant, so not doing any post hoc tests
## save the table
write.table(tests_out, 
            file=paste0(statsFolder, "KS_enhProm_contacts_allGenotypes_allTissues_individualEnhancers.txt"),
            quote=F,
            row.names = F,
            col.names = T)






#################################
########### PLOTTING ############
#################################


# reset the levels to allow nicer plotting
totals$tissue = factor(totals$tissue, levels=c("Un", "Flk1", "CD41"))
totals$genotype = factor(totals$genotype, levels=c("WT", "P1.CTCF.KO", "P2.CTCF.KO"))


## as data are non-normal, plotting using median and IQR!

## E-P interactions split by promoter for Un and Mes

temp = 
  totals %>%
  filter(genotype=="WT" & (tissue=="Un" | tissue=="Flk1"))




## need to identify which points are considered outliers by a standard boxplot and set the limits to that value
# values returned are lower whisker, lower hinge, median, upper hinge, and upper whisker
stat = boxplot.stats(temp$contacts)$stats

lims = c(stat[1],stat[5])

p=
  ggplot(data=temp, aes(x = tissue, y = contacts,fill=tissue))+
  geom_boxplot(size=.1,outlier.shape = NA)+
  coord_cartesian(ylim = lims)+
  facet_wrap(~promName)+
  scale_fill_manual(values=c("#7570b3ff", "#d95f02ff", "#1b9e77ff"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        strip.background =element_rect(fill="white",size=.1),
        axis.ticks = element_line(size=.1),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.1))

ggsave(p, 
       filename = paste0(outFolder, "Un_Mes_WT_median_IQR.pdf"),
       width=2.2,
       height=1.25)


## E-P interactions split by promoter for Mes and CD41

temp = 
  totals %>%
  filter(genotype=="WT" & (tissue=="Flk1" | tissue=="CD41"))

## need to identify which points are considered outliers by a standard boxplot and set the limits to that value
# values returned are lower whisker, lower hinge, median, upper hinge, and upper whisker
stat = boxplot.stats(temp$contacts)$stats

lims = c(stat[1],stat[5])




p=
  ggplot(data=temp, aes(x = tissue, y = contacts,fill=tissue))+
  geom_boxplot(size=.1,outlier.shape = NA)+
  coord_cartesian(ylim = lims)+
  facet_wrap(~promName)+
  scale_fill_manual(values=c("#d95f02ff", "#1b9e77ff"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        strip.background =element_rect(fill="white",size=.1),
        axis.ticks = element_line(size=.1),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.1))


ggsave(p, 
       filename = paste0(outFolder, "Mes_CD41_WT_median_IQR.pdf"),
       width=2.1,
       height=1.25)



## total E-P interactions split by promoter for all tissues and genotypes
temp = 
  totals# %>%
#filter(tissue=="Flk1" | tissue=="CD41")



## need to identify which points are considered outliers by a standard boxplot and set the limits to that value
# values returned are lower whisker, lower hinge, median, upper hinge, and upper whisker
stat = boxplot.stats(temp$contacts)$stats

lims = c(stat[1],stat[5])


p=
  ggplot(data=temp, aes(x = genotype, y = contacts, fill=genotype))+
  geom_boxplot(size=.1,outlier.shape = NA)+
  coord_cartesian(ylim = lims)+
  facet_wrap(promName~tissue, ncol=6)+
  scale_fill_manual(values=c("#666666ff", "#66a61eff", "#e7298aff"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        strip.background =element_rect(fill="white",size=.1),
        axis.ticks = element_line(size=.1),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.1))


ggsave(p, 
       filename = paste0(outFolder, "allTissues_allGenotypes_median_IQR.pdf"),
       width=6,
       height=1.6)



## all enhancers big plot by tissue

p=
  ggplot(data=totals, aes(x = genotype, y = contacts, fill=genotype))+
  stat_summary(geom="bar",
               fun = median)+
  geom_point(alpha=1,size=.3)+
  stat_summary(geom="uperrorbar",
               mapping = aes(x = genotype, y = contacts),
               fun.min = iqr_low,
               fun.max = iqr_hi,
               fun = median,
               width=.3,
               size=.1)+
  facet_wrap(promName~tissue+enh_numName, scales="free", ncol=13)+
  scale_fill_manual(values=c("#666666ff", "#66a61eff", "#e7298aff"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text.x = element_text(size = 6),
        axis.text = element_text(size=6))+
  theme(axis.text.x=element_blank(),
        strip.background =element_rect(fill="white",size=.1),
        axis.ticks = element_line(size=.1),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.1))


ggsave(p, 
       filename = paste0(outFolder, "allTissues_allGenotypes_allEnhancers_median_IQR.pdf"),
       width=13,
       height=17)