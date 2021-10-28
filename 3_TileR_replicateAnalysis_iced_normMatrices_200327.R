### TILED ANALYSIS ###

library(tidyverse)
library(magrittr)

#### taking iced interaction matrices for all non merged individual replicate samples
# and extracting normalised matrices from them to plot


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
#FDRLevel <- 0.1


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



# minumum fragment to view from
minFrag <- 402

# maximum fragment to view from
maxFrag <- 713


###############
### OUTPUTS ###
###############

trackFolder <- paste0(base, "trackFolder_", matrixType, "/")
dir.create(trackFolder, showWarnings = F)

variablesFolder <- paste0(base, "savedVariables_", matrixType, "/")
dir.create(variablesFolder, showWarnings = F)

matrixFolder <- paste0(base, "outMatrixFolder_indiReps_", matrixType, "/")
dir.create(matrixFolder, showWarnings = F)


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

#######################
### DATA WRANGLING  ###
#######################


# assign individual dfs for each sample

for (n in 1:nrow(samples)){
  name=samples[n,1]
  x <- 
    data %>% 
    filter(sampleName==name) %>% 
    mutate(combo=paste0(baitID,"-",preyID)) %>% 
    dplyr::select(-baitID,-preyID)
  assign(paste0("df", n), x)
}


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
  full_join(df9, by = "combo") %>%
  full_join(df10, by = "combo") %>%
  full_join(df11, by = "combo") %>%
  full_join(df12, by = "combo") %>%
  full_join(df13, by = "combo") %>%
  full_join(df14, by = "combo") %>%
  full_join(df15, by = "combo") %>%
  full_join(df16, by = "combo") %>%
  full_join(df17, by = "combo") %>%
  full_join(df18, by = "combo") %>%
  full_join(df19, by = "combo") %>%
  full_join(df20, by = "combo") %>%
  full_join(df21, by = "combo") %>%
  full_join(df22, by = "combo") %>%
  full_join(df23, by = "combo") %>%
  full_join(df24, by = "combo") %>%
  full_join(df25, by = "combo") %>%
  full_join(df26, by = "combo") %>%
  full_join(df27, by = "combo") %>%
  full_join(df28, by = "combo") %>%
  full_join(df29, by = "combo") %>%
  full_join(df30, by = "combo") %>%
  full_join(df31, by = "combo") %>%
  full_join(df32, by = "combo") %>%
  full_join(df33, by = "combo")

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
                      samples[9,1],
                      samples[10,1],
                      samples[11,1],
                      samples[12,1],
                      samples[13,1],
                      samples[14,1],
                      samples[15,1],
                      samples[16,1],
                      samples[17,1],
                      samples[18,1],
                      samples[19,1],
                      samples[20,1],
                      samples[21,1],
                      samples[22,1],
                      samples[23,1],
                      samples[24,1],
                      samples[25,1],
                      samples[26,1],
                      samples[27,1],
                      samples[28,1],
                      samples[29,1],
                      samples[30,1],
                      samples[31,1],
                      samples[32,1],
                      samples[33,1])


## split up the bait and prey IDs
union <- separate(union, `baitID-preyID`, into=c("baitID","preyID"), sep = "-", remove = TRUE)

## UNION IS WIDE FORMAT DATA, POSSIBLY USEFUL so saving
# save this table for use later
saveRDS(union, paste0(variablesFolder, "/Tiled_indiReps_", matrixType, ".rds"))
write.table(union, file=paste0(variablesFolder, "/Tiled_indiReps_", matrixType, ".txt"),
            col.names = T, 
            row.names = F,
            quote=F,
            sep="\t")

## to load it back
#union <- readRDS(paste0(variablesFolder, "/Tiled_indiReps_", matrixType, ".rds"))


########################
### NORMALISE SIGNAL ###
########################

# make all columns of union numeric

for (i in 1:length(union)){
  union[,i] <- as.numeric(union[,i])
}


# before plotting matrices, normalising between samples to the mean number of reporter contacts
# so half samples will be scaled up, and half will be scaled down

# mean is 3 990 509 (from Numbers.txt)

data.matrix <- data.matrix(union) 
column.totals <- apply(union, 2, sum)
data.norm <- t(t(data.matrix) * 3990509 / column.totals)
norm.totals <- apply(data.norm, 2, sum)
column.totals
norm.totals

# combine with union
normed <- cbind(union[,1:2],
                data.norm[,3:ncol(data.norm)])


# rename columns
# this finds combined and replaces with combined_normed
colnames(normed) <- gsub("(combined)","\\comb_normed",colnames(normed))


## UNION IS WIDE FORMAT DATA, POSSIBLY USEFUL so saving
# save this table for use later
saveRDS(normed, paste0(variablesFolder, "/Tiled_indiReps_normed_", matrixType, ".rds"))
write.table(normed, file=paste0(variablesFolder, "/Tiled_indiReps_normed_", matrixType, ".txt"),
            col.names = T, 
            row.names = F,
            quote=F,
            sep="\t")



# export the normalised matrices
# remove any values with a zero to decrease file size
for (i in 1:length(normed)){
  a = i+2
  
  if (a > length(normed)){
    break
  } else {
   
    export = 
    normed %>%
      dplyr::select(1,2,a) %>%
      filter(.[,3]>0)
    
    exp_file = paste0("Tiled-C_matNorm_", colnames(export)[3], ".matrix")
    
    write.table(export,
                file=paste0(matrixFolder,exp_file),
                quote=F,
                row.names=F,
                col.names=F)
  }
}



###########################################################



########################################
### FILTER NORMALISED DATA ON WINDOW ###
########################################
# as in a smalling number of bins to view


# make into integers to allow filtering
fragData$baitID <- as.integer(fragData$baitID)
fragData$preyID <- as.integer(fragData$preyID)
normed$baitID <- as.integer(normed$baitID)
normed$preyID <- as.integer(normed$preyID)


cat("Now selecting just a window of the data...", "\n")
  
# filter the union or normed table on just the desired minimum and maximum fragments
filtData <- 
  normed %>%
  filter(baitID > minFrag &
         baitID < maxFrag &
         preyID > minFrag &
         preyID < maxFrag)

# for plotting in python, I need the bins numbered 0 onwards!
# to get this just subtract the minFrag from each
# these are then relative to minFrag=0
# -1 to correct for 0 indexing
filtData %<>%
mutate(relBaitID=baitID-minFrag-1,
       relPreyID=preyID-minFrag-1)

### finding out the viewing window to use on UCSC
UCSCview <- 
  fragData %>%
  filter(baitID > minFrag &
           baitID < maxFrag &
           preyID > minFrag &
           preyID < maxFrag)

chr <- UCSCview$prey_chr[1]
UCSCmin <- min(UCSCview$prey_start)
UCSCmax <- max(UCSCview$prey_end)
binNum <- nrow(UCSCview)


# export the normalised matrices
# remove any values with a zero to decrease file size
for (i in 1:length(normed)){
  a = i+2
  
  if (a > length(normed)){
    break
  } else {
    
    export = 
      filtData %>%
      dplyr::select(starts_with("rel"),a) %>%
      filter(.[,3]>0)
    
    exp_file = paste0("Tiled-C_matNorm_window_", chr, "-", UCSCmin,"-", UCSCmax,"_", binNum, "_bins_", colnames(export)[3], ".matrix")
    
    write.table(export,
                file=paste0(matrixFolder,exp_file),
                quote=F,
                row.names=F,
                col.names=F)
  }
}

  
cat("Export of window filtered normalised matrix for viewing window", minFrag, "to", maxFrag, "comlete", "\n")

