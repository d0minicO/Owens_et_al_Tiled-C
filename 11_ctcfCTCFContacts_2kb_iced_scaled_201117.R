### TILED ANALYSIS ###

library(tidyverse)
library(tools)
library(cowplot)
library(magrittr)
library(RColorBrewer)
library(broom)
library(plotrix)

library(FSA)
library(ggpubr)

#### taking iced and scaled interaction matrices calculated by Marieke's script 
# from individual samples

# calculating CTCF CTCF contacts just for TAD boundaries in different groups (with replicate info as SD)


# created 18 11 2020
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


## promData
# promoter data including location of bin containing CTCF site deletions
# need this for excluding it later
# format is following cols in tab separated fle
# chr start stop  promName  promBin CTCF_dels_bin
promFile <- paste0(base, "promData.txt")


## ctcfData
# had bed file of ctcf peak cordinates and using bedtools 
# intersected with fragData to get overlapping bins
ctcfFile <- paste0(base, "ctcfData_bins_clip.txt")


###############
### OUTPUTS ###
###############

outFolder <- paste0(base, "ctcfCTCF_contacts_", matrixType, "/")
dir.create(outFolder, showWarnings = F)

variablesFolder <- paste0(base, "savedVariables_", matrixType, "/")
dir.create(variablesFolder, showWarnings = F)



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


###########################
### LOAD THE CTCF sites ###
###########################

## load the ctcf data
ctcfData=read_tsv(ctcfFile,col_names=c("chr", "start", "stop", "ctcfID", "chr_ctcf", "start_ctcf", "stop_ctcf", "ctcfName"))


# filter to make the rows unique
ctcfData_new <-
  ctcfData %>%
  distinct(ctcfID, .keep_all = T) %>%
  arrange(ctcfID)

# add a better just numerical ctcf name to use for plotting later
Runx1_ATG=92825890

ctcfData_new %<>%
  mutate(ctcf_mid=round((start_ctcf+stop_ctcf)/2)) %>%
  mutate(ctcf_numName=factor(round(((Runx1_ATG-ctcf_mid)/1e3),digits=1)))

# output the ctcf data table
write.table(ctcfData_new,
            file=paste0(outFolder,"CTCFData_all_clip.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            na = "",
            sep="\t")




### now actually manually clipped the CTCF sites just to the boundary ones
## four outer most CTCF sites
ctcfData_new = 
  read_tsv(file=paste0(outFolder,"CTCFData_all_clip_boundaries.txt"),
           col_names=c("chr", "start", "stop", "ctcfID", "chr_ctcf", "start_ctcf", "stop_ctcf", "ctcfName", "ctcfMid", "ctcf_numName"))





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



####################################
#### GATHER CTCF CTCF CONTACTS #####
####################################
totals = data.frame()

for (i in 1:nrow(ctcfData_new)){ # for ctcf loop
  
  ctcfName=ctcfData_new$ctcf_numName[i]
  ctcfID=ctcfData_new$ctcfID[i]
  
  
  # filter the matrix on just this ctcf
  ctcf_filt <-
    d %>% 
    filter(baitID==ctcfID | preyID==ctcfID) %>%
    mutate(ctcfName=ctcfName)
  
  
  # filter this matrix on all other CTCF sites
  ctcf_filt2 =
    ctcf_filt %>%
    filter(baitID %in% ctcfData_new$ctcfID & preyID %in% ctcfData_new$ctcfID)
  
  
  # combine to totals df
  totals=rbind.data.frame(totals,ctcf_filt2)
  
} # end of for prom loop



## remove the weird instances where I am quantifying bins interacting with themselves!
totals %<>%
  filter(baitID!=preyID)



## subset on just the long range boundary interactions (ie the top of the TAD triangle)
totals %<>%
  mutate(binDiff=preyID-baitID) %>%
  filter(binDiff>100)



# get total CTCF contacts in each group
# start by getting it for each sample
# and then split up the sample into group
sum_CTCF_contacts <-
  totals %>%
  group_by(sample,ctcfName) %>%
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
          main = "Density plot of per CTCF site CTCF interactions",
          xlab = "per per CTCF site CTCF interactions")


ggqqplot(totals$contacts)

## both look highly non-normal

## doing a shapiro-wilks test of normality
# had to sample as max for func is 5k
shapiro.test(totals$contacts)
#W = 0.89157, p-value < 2.2e-16

## data is non-normal and so non-parametric Kruskall Wallis tests should be used



### TEST 1 ###

## WILD TYPE SAMPLES ACROSS DIFFERENTIATION

## I know beforehadn that I am only interested in WT samples now,
# so it is okay to do stat test on this smaller subset of the data
# adjusted p values will be given in the Dunn multiple comparisons test

## WT only
temp=
  totals %>% 
  filter(genotype=="WT")

t1 = kruskal.test(contacts ~ tissue, data = temp)
print(t1)

#Dunn's Kruskal-Wallis post-hoc test
posthoc1<-dunnTest(contacts ~ tissue, data=temp, method=meth)
ph1 = 
  print(posthoc1) %>%
  mutate(Dunn_sig=if_else(P.adj<0.05, "SIG*", "not"),
         testGroup=paste0("WT"),
         KW_pval=tidy(t1)$p.value,
         KW_signif=if_else(tidy(t1)$p.value<0.05, "SIG*", "not"))


## save the table
write.table(ph1, 
            file=paste0(outFolder, "KW_Dunn_CTCF-TAD-boundaries_WT_allTissues.txt"),
            quote=F,
            row.names = F,
            col.names = T)




#################################
########### PLOTTING ############
#################################


# reset the levels to allow nicer plotting
totals$tissue = factor(totals$tissue, levels=c("Un", "Flk1", "CD41"))
totals$genotype = factor(totals$genotype, levels=c("WT", "P1.CTCF.KO", "P2.CTCF.KO"))


### PLOT 1 ###

# WT just Un and Mes MAIN TAD ONLY

temp = 
  totals %>%
  filter(genotype=="WT" & (tissue =="Un" | tissue=="Flk1"))

# check the filtering work a I had some issues
levels(factor(temp$tissue))
levels(factor(temp$genotype))


#check what the median values will be
temp %>%
  group_by(tissue) %>%
  summarise(median(contacts))


## need to identify which points are considered outliers by a standard boxplot and set the limits to that value
# values returned are lower whisker, lower hinge, median, upper hinge, and upper whisker
stat = boxplot.stats(temp$contacts)$stats

lims = c(stat[1],stat[5]*1.1)




p=
  ggplot(data=temp, aes(x = tissue, y = contacts, fill=tissue))+
  geom_boxplot(size=.1,outlier.shape = NA)+
  stat_summary(fun.data = n_fun, geom = "text",color="black",size=1.5,vjust=-1,aes(y=stat[5]))+
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
       filename = paste0(outFolder, "Un_Mes_WT_CTCF-TAD-boundaries_median_IQR.pdf"),
       width=1.7,
       height=1)




### PLOT 2 ###

# WT just MES and HPC MAIN TAD ONLY

temp = 
  totals %>%
  filter(genotype=="WT" & (tissue =="Flk1" | tissue=="CD41"))

# check the filtering work a I had some issues
levels(factor(temp$tissue))
levels(factor(temp$genotype))


#check what the median values will be
temp %>%
  group_by(tissue) %>%
  summarise(median(contacts))



## need to identify which points are considered outliers by a standard boxplot and set the limits to that value
# values returned are lower whisker, lower hinge, median, upper hinge, and upper whisker
stat = boxplot.stats(temp$contacts)$stats

lims = c(stat[1],stat[5]*1.1)




p=
  ggplot(data=temp, aes(x = tissue, y = contacts, fill=tissue))+
  geom_boxplot(size=.1,outlier.shape = NA)+
  stat_summary(fun.data = n_fun, geom = "text",color="black",size=1.5,vjust=-1,aes(y=stat[5]))+
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
       filename = paste0(outFolder, "Mes_HPC_WT_CTCF-TAD-boundaries_median_IQR.pdf"),
       width=1.8,
       height=1)


