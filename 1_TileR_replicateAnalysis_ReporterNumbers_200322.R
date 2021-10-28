### TILED ANALYSIS ###

library(tidyverse)
library(magrittr)
library(ggplot2)
library(RColorBrewer)

#### taking iced interaction matrices for all samples and calculating reporter numbers from them


#################
### FUNCTIONS ###
#################



##############
### INPUTS ###
##############


## the main directory
base = "C:/Users/Dominic/Desktop/Work/Paper/Bioinformatics/Tile-C/CTCF-KO/"



####  relative to base folder, or could be specified elsewhere

## the bait file
baitFile <- paste0(base, "baitsRightIDs.txt")

## the fragments genome file
fragFile <- paste0(base, "fragData_2kb.txt")

## the data folder
dataFolder <- paste0(base, "iced_matrix/")


###############
### OUTPUTS ###
###############

## DESeq directories
QC_dir <-  paste0(base, "QC/")
dir.create(QC_dir, showWarnings=F)



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
# clone names (reps effectively)
clones <- levels(factor(word(samples$sample, 3, sep="_")))



###########################
#### DATA GATHER LOOP #####
###########################


#initialise empty data table
data <- NULL

# loop over directories and files to gather lots of files

for (i in 1:nrow(samples)){ # for samples loop
  sampleFile = samples[i,2]
  
  # get the variable for which sample
  thisSample <- samples[i,1]
  
  # load the data and add a column for the sample
  interactions <- 
    data.table::fread(sampleFile, fill = TRUE) %>%
    mutate(sample=thisSample)
  
  # combine all in one df
  data <- rbind(data,interactions)
} # for samples loop


colnames(data) <- c("baitID", "preyID", "reads", "sampleName")

cat("Total Usable Reporters:", sum(data$reads), "\n")

# Total Usable Reporters: 131686783


# total reporters in each sample
reportersBySample <-
  data %>%
  group_by(sampleName) %>%
  summarise(total=sum(reads))


# mean and sd of reporters in each sample
reportersBySample %>%
  summarise(mean=mean(total),
            sd=sd(total)) 


# organise into IDs and sort based on it
reportersBySample %<>%
  mutate(sampleName2=sampleName) %>%
  separate(sampleName2,into=c(NA,NA,NA,NA,"ID"),sep="_") %>%
  separate(ID,into=c(NA,"ID"),sep="ID") %>%
  mutate(ID=as.numeric(ID)) %>%
  arrange(ID)

write_tsv(reportersBySample,path = paste0(QC_dir, "Numbers.txt"), append = T, col_names = T)


# give the 3C libs genotype and cellType variables to help plotting
reportersBySample %<>%
  mutate(sampleName2=sampleName) %>%
  separate(sampleName2,into=c("cellType","genotype"),sep="_")


# plot just the reporter numbers
p<-
  ggplot(reportersBySample, aes(x=sampleName,y=total))+
  geom_bar(stat="identity", aes(fill=cellType,colour=genotype),size=.5,width=.5)+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Blues")+
  xlab(label="Tiled-C Sample") +
  ylab(label="Total reporter counts")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))


ggsave(p,
       filename=paste0(QC_dir,"/Total_IndiRep_reporter_counts.pdf"),
       width=7,
       height=5)
  
  # will later correlate to cell numbers etc