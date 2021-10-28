### TILED ANALYSIS ###

library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(gplot)
library(ggplot2)
library(stringr)

#### This script is to calculate interaction frequencies per 2kb window and plots them 
## all done on RAW CONTACT MATRICES

# working 04 04 2020



###############
### OPTIONS ###
###############

options(max.print=50)



#################
### FUNCTIONS ###
#################


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

# the significance level to use for DESeq
FDRLevel <- 0.1


# are we using iced or raw matrices?
matrixType = "raw_merged" # "raw"


####  relative to base folder, or could be specified elsewhere

## the bait file

###### NOTE HERE WE ARE USING ALL THE BAITS IN THE ENTIRE TILED REGION ######
baitFile <- paste0(base, "fragData_2kb.txt")

## the fragments genome file
fragFile <- paste0(base, "fragData_2kb.txt")

## the data folder
if (matrixType=="iced"){
  dataFolder <- paste0(base, "iced_matrix/")
} else if (matrixType=="raw") {
  dataFolder <- paste0(base, "matrix/")
} else if (matrixType=="raw_merged") {
  
  dataFolder <- paste0(base, "merged_matrix/")
} else {
  cat("Warning: matrixType must be set to either \"iced\" or \"raw\"")
}



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

baits %<>%
  mutate(baitName=V4)

colnames(baits) <- c("bait_chr", "bait_start", "bait_end", "baitID","baitName")
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


# get the matrix location (sampleFIle) and add to the table
samples %<>%
  mutate(sampleFiles=paste0(dataFolder,sample,"/raw/2000/tiled_",sample,"_2000.matrix"))

# get the matrix location (sampleFIle) and add to the table
# along with
# tissue types, genotypes, and clone names
samples %<>%
  mutate(sampleFiles=paste0(dataFolder,sample,"/raw/2000/tiled_",sample,"_2000.matrix"),
         sample2=sample) %>%
  separate(sample2,into=c("tissue","genotype","clone"),sep="_")



# get the tissue types
tissues <- levels(factor(word(samples$sample, 1, sep="_")))
# geno types
genotypes <- levels(factor(word(samples$sample, 2, sep="_")))
# exchange the - for . to prevent problems with column naming
genotypes <- gsub("-", ".", genotypes)

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


#######################
### DATA WRANGLING  ###
#######################

# now do a full join by bait and prey IDs combo 
# so that any sample with a read for that bait-prey combo 
# will be assigned to it in each row
union <- 
  df1 %>% 
  full_join(df2, by = "combo") %>%
  full_join(df3, by = "combo") %>%
  full_join(df4, by = "combo") %>%
  full_join(df5, by = "combo") %>%
  full_join(df6, by = "combo") %>%
  full_join(df7, by = "combo") %>%
  full_join(df8, by = "combo") %>%
  full_join(df9, by = "combo")

# keep only the reads
union <-
  union %>%
  dplyr::select(-starts_with("sampleName")) %>%
  dplyr::select(combo,starts_with("reads"))

# convert NAs to 0
union[is.na(union)] <- 0

# make sample names not a factor
samples$sample <- as.character(samples$sample)

# give proper names to the read columns
colnames(union) <- c("baitID-preyID",
                     samples[1,1],
                     samples[2,1],
                     samples[3,1],
                     samples[4,1],
                     samples[5,1],
                     samples[6,1],
                     samples[7,1],
                     samples[8,1],
                     samples[9,1])

## split up the bait and prey IDs
union <- separate(union, `baitID-preyID`, into=c("baitID","preyID"), sep = "-", remove = TRUE)

# make into integers to allow joining
union$baitID <- as.integer(union$baitID)
union$preyID <- as.integer(union$preyID)
fragData$baitID <- as.integer(fragData$baitID)
fragData$preyID <- as.integer(fragData$preyID)



## alternatively moving on for plotting and normalisation  #################
# join on bait and prey IDs
# only keeping the frags in format chr:start-stop
union <- 
  union %>% 
  full_join(fragData, by = "baitID") %>%
  dplyr::rename(preyID=preyID.x) %>%
  full_join(fragData, by = "preyID") %>%
  dplyr::select(1,12:14,2,17:19,3:11) %>%
  mutate(#"comboID"=paste0(baitID.x,"-",preyID),
         "bait_frag"=paste0(prey_chr.x,":",prey_start.x,"-",prey_end.x),
         "prey_frag"=paste0(prey_chr.y,":",prey_start.y,"-",prey_end.y)) %>%
  dplyr::rename(baitID=baitID.x) %>%
  dplyr::select(-matches("\\.x"),-matches("\\.y")) %>%
  dplyr::select(baitID,preyID,bait_frag,prey_frag,3:11)

#union2=union

##########################
### BAIT SPECIFIC PART ###
##########################

# calcuating RAW interaction count per bait per replicate


# initialise empty variables for calculating interaction number per bait
# for raw and scaled normalised to read count
intNum <- data.frame()
normIntNum <- data.frame()

for (i in 1:nrow(baits)){ # working on each bait at a time
  
  # get the bait_frag, baitName, and baitChrom we are working on
  baitUsed = baits$bait_frag[i]
  baitName = baits$baitName[i]
  baitNum = baits$baitID[i]
  #baitChrom = 
  #  baits[i,] %>%
  #  separate(bait_frag, into=c("baitChrom", NA), sep=":", remove=T) %>%
  #  dplyr::select(baitChrom) %>%
  #  as.character()
  
  cat("Now analysing bait", baits$baitName[i], "Tile-C version", "\n")
  
  # filter the union table on just this bait
  # note I am actually not really using the concept of bait and prey
  # just looking for any fragment that interacts with this viewpoint!!
  # grep on the combo of "baitID-preyID" wasn't working as grep(187) found 1871 etc
  # also need to get a fragID column to plot the data over
  # also excluding the baited bin
  
  countData <- 
    union %>%
    filter(baitNum==baitID | baitNum==preyID) %>%
    mutate(fragID=((baitID+preyID)-baitNum)) %>%
    filter(fragID!=baitNum) %>%
    dplyr::select(-baitID,-preyID,-bait_frag,-prey_frag) %>%
    tibble::column_to_rownames(var = "fragID")
    

  #####################
  ### NORMALISATION ###
  #####################
  
  # make all columns of countData numeric to allow summing
  
  for (i in 1:length(countData)){
    countData[,i] <- as.numeric(countData[,i])
  }
  
  # exact normalisation method from CC stats pipeline
  data.matrix <- data.matrix(countData) 
  column.totals <- apply(countData, 2, sum)
  data.norm <- t(t(data.matrix) * 5e3 / column.totals) 
  column.totals.norm <- apply(data.norm, 2, sum)
  
  # add the column totals to intNum
  
  intNum <- rbind.data.frame(intNum,column.totals)
  normIntNum <- rbind.data.frame(normIntNum,column.totals.norm)

  

}




#######################################
### EXPLORING INTERACTIONS PER BAIT ###
#######################################

# rename the outputted tables columns
colnames(intNum) <- samples$sample
colnames(normIntNum) <- samples$sample

# export wide format table
write.table(intNum, 
            file=paste0(QC_dir,"Tiled_2kb_raw_merged_intNumPer2kbFragPerSamples_wide.txt"),
            quote=F,
            col.names=T,
            row.names=F)



# make into long format
long <- gather(data = intNum,
               key = "sample",
               value = "contacts",
               1:9,
               factor_key = T)


# get the cellType_genoType group
long <-
  long %>%
  mutate(sample2=sample) %>%
  separate(sample2, into=c("cellType", "genoType","combo"),sep="_") %>%
  mutate(genoType=gsub("-",".",genoType)) %>%
  mutate(group=paste0(cellType,"_",genoType))

# convert to factor
long$group <- factor(long$group)

# get the numbers in each library
meanCont <- long %>%
  group_by(sample) %>%
  summarise(mean(contacts))

sdCont <- long %>%
  group_by(sample) %>%
  summarise(sd(contacts))

minCont <- long %>%
  group_by(sample) %>%
  summarise(min(contacts))

maxCont <- long %>%
  group_by(sample) %>%
  summarise(max(contacts))


# get the numbers in each merged group (cellType_genoType)
sumContCell <- long %>%
  group_by(group) %>%
  summarise(sum(contacts))

meanContCell <- long %>%
  group_by(group) %>%
  summarise(mean(contacts))

sdContCell <- long %>%
  group_by(group) %>%
  summarise(sd(contacts))


# get the numbers across all libraries
meanContAll <- long %>%
  summarise(mean(contacts))

sdContAll <- long %>%
  summarise(sd(contacts))


## quick plotting

box <- 
  ggplot(long, aes(x=sample, y=contacts, fill=cellType,colour=genoType)) +
    geom_boxplot() +
    scale_color_brewer(palette="Dark2")+
    scale_fill_brewer(palette="Blues")+
    xlab(label="Tiled-C Sample") +
    ylab(label="Reporter counts per bin")+
  ggtitle("Reporter counts per bin by Tiled-C library")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=8))


ggsave(plot=box, 
       filename=paste0(QC_dir,"Average_Tiled-C_repCountsPerBin_2kb_raw_merged.pdf"),
       device="pdf",
       width=4.5,
       height=4)