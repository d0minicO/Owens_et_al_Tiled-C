### TILED ANALYSIS ###

library(tidyverse)
library(tools)
library(cowplot)
library(magrittr)
library(RColorBrewer)

#### taking iced and scaled interaction matrices calculated by Marieke's script 
# from individual samples

# calculating enhancer promoter contacts in different groups (with replicate info as SD)


# working 23 04 2020


###############
### OPTIONS ###
###############

options(max.print=50)


#################
### FUNCTIONS ###
#################

head_tail <- function(x,description){
  # this function takes a data frame (x) and adds a header containing the "description" string
  # header is added to first row of x
  # tail is added of two empty rows
  
  rbind(c(description, rep("",ncol(x)-1)),
        x,
        rep("",ncol(x)),
        rep("",ncol(x)))
  
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



###############
### OUTPUTS ###
###############

outFolder <- paste0(base, "enhProm_contacts_", matrixType, "/")
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


## load the promData
promData=read_tsv(promFile,col_names = c("chr", "start", "stop", "promID", "promName", "delID"))

# get the IDs of the bins with actual CTCF site deletions in to exclude from the counting
# as these would not allow for a fair comparison between genotypes (within genotype would be fine though)
delIDs=promData$delID


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
# and check it all makes sense

# iced individual replicate data was previously anlaysed in 3_TileR_replicateAnalysis_iced_normMatrices_200327
#normed <- readRDS(paste0(variablesFolder, "/Tiled_indiReps_normed_", matrixType, ".rds"))

# check that normalisation to mean contacts in each sample worked (3 990 509 from Numbers.txt)
#apply(normed, 2, sum)

# update the col names
# exchange the - for . to prevent problems with column naming
#colnames(normed) <- gsub("-", ".", colnames(normed))

# get into long format
# split up sample to get cellType, genoType, Clone, XP, and ID cols
#d <- normed %>%
#  gather(sample, contacts, -c(baitID, preyID)) %>%
#  mutate(sample2=sample) %>%
#  separate(sample2,into=c("tissue","genotype", "clone", "XP", "ID"),sep="_") %>%
#  mutate(group=paste0(tissue,"_",genotype))
  

# check that the longform and grouping still works
#d %>%
#  group_by(sample) %>%
#  summarise(sum=sum(contacts))

#CD41_P1-CTCF-KO_E7-1_XP1_ID5   3990509 
#2 CD41_P1-CTCF-KO_E7-1_XP4_ID26  3990509.
#3 CD41_P1-CTCF-KO_E7-2_XP2_ID14  3990509 
#etc

#d %>%
#  group_by(group) %>%
#  summarise(sum=sum(contacts))

#1 CD41_P1-CTCF-KO 15962036
#2 CD41_P2-CTCF-KO 15962036
#3 CD41_WT         15962036
#4 Flk1_P1-CTCF-KO 15962036
#5 Flk1_P2-CTCF-KO 15962036
#6 Flk1_WT         15962036
#7 Un_P1-CTCF-KO   11971527
#8 Un_P2-CTCF-KO   11971527
#9 Un_WT           11971527


# save this data
#saveRDS(d, paste0(variablesFolder, "/Tiled_indiReps_normed_", matrixType, "_longFormat.rds"))

# load previously processed data instead of processing again as it takes a while!
d <- readRDS(paste0(variablesFolder, "/Tiled_indiReps_normed_", matrixType, "_longFormat.rds"))


#########################################
#### GATHER TOTAL PROMOTER CONTACTS #####
#########################################
totals = data.frame()

for (i in 1:nrow(promData)){ # for prom loop
  
  promName=promData$promName[i]
  promID=promData$promID[i]

  
  # filter the matrix on just this prom
  prom_filt <-
    d %>% 
    filter(baitID==promID | preyID==promID)
  
  # filter out all the excluded delIDs
  for (delID in delIDs){
    prom_filt %<>% 
      filter(!baitID %in% delID &
               !preyID %in% delID)
  }

  
  # get total promoter contacts in each group
  # start by getting it for each sample
  # and then split up the sample into group
  sum_promContacts <-
    prom_filt %>%
    group_by(sample) %>%
    summarise(sum=sum(contacts)) %>%
    separate(sample,into=c("tissue","genotype", "clone", "XP", "ID"),sep="_") %>%
    mutate(group=paste0(tissue,"_",genotype),
           promName=promName) %>%
    ungroup()
  
  
  # add to output df
  totals=rbind.data.frame(totals,sum_promContacts)
  
} # end of for prom loop

###########################################
#### ANOVA ON TOTAL PROMOTER CONTACTS #####
###########################################

# do stats (three way anova)
res.aov <- aov(sum ~ tissue*genotype*promName, data = totals)
# Summary of the analysis
summary(res.aov)

# posthoc test
# rest was useful when had interaction effects to explore
postHoc <-
  data.frame(TukeyHSD(res.aov)$promName) %>%
  rownames_to_column("comparison")# %>%
  #filter(p.adj<0.1) %>%
  #dplyr::select(comparison,diff,p.adj) %>%
  #mutate(diff=signif(diff,3),
  #       p.adj=signif(p.adj,3)) %>%
  #tidyr::separate(comparison,into=c("group1", "group2"),sep="-") %>%
  #separate(group1, into=c("cellType1","genoType1"),sep="_") %>%
  #separate(group2, into=c("cellType2","genoType2"),sep="_") %>%
  #separate(genoType1, into=c("genoType1", "prom1"),sep=":") %>%
  #separate(genoType2, into=c("genoType2", "prom2"),sep=":")


# filter for the different significant within group comparisons
proms_comp =
  postHoc# %>%
#  filter(cellType1==cellType2 & genoType1==genoType2)

# no sigDifferences between genotypes
#genoType_comp =
#  postHoc %>%
#  filter(cellType1==cellType2 & prom1==prom2)

#cellType_comp =
#  postHoc %>%
#  filter(genoType1==genoType2 & prom1==prom2)


## add header and tail
proms_comp = head_tail(proms_comp,"Comparing promoters")
#genoType_comp = head_tail(genoType_comp,"Comparing genotypes")
#cellType_comp = head_tail(cellType_comp,"Comparing cell types")

# combine altogether and save
postHoc_out <- rbind.data.frame(proms_comp)
write.table(postHoc_out,
            file=paste0(outFolder,"Three-Way-ANOVA_Tukeys-PostHoc_allComparisons.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            sep="\t")


#######################################
#### PLOT TOTAL PROMOTER CONTACTS #####
#######################################

# now calculate mean and sd of the sum of contacts for each group
output <- 
  totals %>%
  group_by(tissue,genotype,promName) %>%
  summarise(mean=mean(sum),
            sd=sd(sum))


# plot these by cell type, genotype, and promoter
p1 <-
  ggplot(output, aes(x=genotype,y=mean))+
  geom_bar(stat = "identity",position=position_dodge(), aes(fill=genotype))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),position=position_dodge(.9),width=.2)+
  scale_fill_brewer(palette="Accent")+
  facet_grid(~tissue+promName)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))


p2 <-
  ggplot(output, aes(x=tissue,y=mean))+
  geom_bar(stat = "identity",position=position_dodge(), aes(fill=tissue))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),position=position_dodge(.9),width=.2)+
  scale_fill_brewer(palette="Accent")+
  facet_grid(~genotype+promName)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

p3 <-
  ggplot(output, aes(x=promName,y=mean))+
  geom_bar(stat = "identity",position=position_dodge(), aes(fill=promName))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),position=position_dodge(.9),width=.2)+
  scale_fill_brewer(palette="Accent")+
  facet_grid(~genotype+tissue)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

# save plot
combo = plot_grid(p1, p2, p3, align="v",ncol=1,axis="rlbt")
ggsave(plot=combo,
       filename=paste0(outFolder, "totalPromContacts.pdf"),
       width=9,
       height = 7)



# plot just the different promoters across cell type/genotypes
just_proms <-
  totals  %>%
  ungroup() %>%
  group_by(promName) %>%
  summarise(mean=mean(sum),
            sd=sd(sum))


p4 <-
  ggplot(just_proms, aes(x=promName,y=mean))+
  geom_bar(stat = "identity",position=position_dodge(), aes(fill=promName))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),position=position_dodge(.9),width=.2)+
  scale_fill_brewer(palette="Accent")+
  #facet_grid(~genotype)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))


ggsave(plot=p4,
       filename=paste0(outFolder, "totalPromContacts_justProms.pdf"),
       width=3,
       height = 2.5)
