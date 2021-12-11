### TILED ANALYSIS ###

library(tidyverse)
library(magrittr)

#### taking iced interaction matrices for merged samples
# Normalising to mean reporter counts across all samples, performing subtraction between the cell types
# subsetting the interactions for a given "zoom" (window)
# also filtering on DESeq2 differentially interacting fragments calculated in 
# TileR_replicateAnalysis_DESeq_diffIntFrags_nonIced_matrices_v2_200324

# CTCF-KO version 27nd March 2020

###############
### OPTIONS ###
###############

options(max.print=25)

#################
### FUNCTIONS ###
#################


##############
### INPUTS ###
##############


## the main directory
base = "C:/Users/Dominic/Desktop/Work/Paper/Bioinformatics/Tile-C/CTCF-KO/"

#### relative to base folder, or could be specified elsewhere

## the data folder
## MERGED (n=3 or 4) ICED MATRICES
dataFolder <- paste0(base, "merged_iced_matrix/")

# are we using iced or raw matrices?
matrixType = "iced" # "raw"

## the bait file
baitFile <- paste0(base, "baitsRightIDs.txt")

## the fragments genome file
fragFile <- paste0(base, "fragData_2kb.txt")



# minumum fragment to view from
# this is for standard window over RUnx1 gene
# that was used in my thesis
# to make 416B and E14 windowed matrix plots
#minFrag <- 402

# maximum fragment to view from
#maxFrag <- 713


## EOMES window
minFrag <- 499
maxFrag <- 735



# DESeq2 output folder that contains the relevant diffInt frags to subset matrices on
FDRLevel <- 0.1
DESeqDir <-  paste0(base, "DESeq_FDR_", FDRLevel, "_raw/")
outDirData = paste0(DESeqDir, "data_raw/")


###############
### OUTPUTS ###
###############

matrixFolder <- paste0(base, "outMatrixFolder_merged_iced/")
dir.create(matrixFolder, showWarnings = F)

variablesFolder <- paste0(base, "savedVariables_", matrixType, "/")
dir.create(variablesFolder, showWarnings = F)


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

## UNION IS WIDE FORMAT DATA, POSSIBLY USEFUL so saving
# save this table for use later
saveRDS(union, paste0(variablesFolder, "/Tiled_merged_", matrixType, ".rds"))
write.table(union, file=paste0(variablesFolder, "/Tiled_merged_", matrixType, ".txt"),
            col.names = T, 
            row.names = F,
            quote=F,
            sep="\t")

## to load it back
#union <- readRDS(paste0(variablesFolder, "/Tiled_merged_", matrixType, ".rds"))


########################
### NORMALISE SIGNAL ###
########################

# make all columns of union numeric

for (i in 1:length(union)){
  union[,i] <- as.numeric(union[,i])
}


# before subtracting or plotting together, normalising between samples to the mean number of reporter contacts
# so half samples will be scaled up, and half will be scaled down

# mean is 14  631 865 (from Numbers_mergedLibs.txt)

data.matrix <- data.matrix(union) 
column.totals <- apply(union, 2, sum)
data.norm <- t(t(data.matrix) * 14631865 / column.totals)
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
saveRDS(normed, paste0(variablesFolder, "/Tiled_merged_normed_", matrixType, ".rds"))
write.table(normed, file=paste0(variablesFolder, "/Tiled_merged_", matrixType, ".txt"),
            col.names = T, 
            row.names = F,
            quote=F,
            sep="\t")


# load it back
normed <- readRDS(paste0(variablesFolder, "/Tiled_merged_normed_", matrixType, ".rds"))



# export the normalised matrices
for (i in 1:length(normed)){
  a = i+2
  
  if (a > length(normed)){
    break
  } else {
   
    export = 
    normed %>%
      dplyr::select(1,2,a) %>%
      filter(.[,3]!=0)
    
    exp_file = paste0("Tiled-C_matNorm_", colnames(export)[3], ".matrix")
    
    write.table(export,
                file=paste0(matrixFolder,exp_file),
                quote=F,
                row.names=F,
                col.names=F)
  }
}



###########################################################


###########################
### COMPARING GENOTYPES ###
###########################

# do all possible subtracts between GENOTYPES on the normalised data
# ie within cell types


for (tissue in tissues){
  
  subset =
    normed %>%
      select(1,2,starts_with(tissue))
  
  sbtrct <-
    normed %>%
    mutate(subAminC=subset[,3]-subset[,5],
           subBminC=subset[,4]-subset[,5],
           subAminB=subset[,3]-subset[,4],
           subCminA=subset[,5]-subset[,3],
           subCminB=subset[,5]-subset[,4]) %>%
    select(1,2,starts_with("sub"))
  
  colnames(sbtrct)[3:7] = c(paste(tissue,genotypes[1],"minus",genotypes[3],sep="_"),
                            paste(tissue,genotypes[2],"minus",genotypes[3],sep="_"),
                            paste(tissue,genotypes[1],"minus",genotypes[2],sep="_"),
                            paste(tissue,genotypes[3],"minus",genotypes[1],sep="_"),
                            paste(tissue,genotypes[3],"minus",genotypes[2],sep="_"))
  

  # export the subtracted matrices
  for (i in 1:length(sbtrct)){
    a = i+2
    
    if (a > length(sbtrct)){
      break
    } else {
      
      export = 
        sbtrct %>%
        dplyr::select(1,2,a) %>%
        filter(.[,3]!=0)
      
      exp_file = paste0("Tiled-C_matNorm_subtract_", colnames(export)[3], ".matrix")
      
      write.table(export,
                  file=paste0(matrixFolder,exp_file),
                  quote=F,
                  row.names=F,
                  col.names=F)
    }
  }
}



#########################
### COMPARING TISSUES ###
#########################
# now subtracts will be made within GENOTYPES but accross cell types


for (genotype in genotypes){
  
  subset =
    normed %>%
    select(1,2,matches(genotype))
  
  sbtrct <-
    normed %>%
    mutate(subAminC=subset[,3]-subset[,5],
           subBminC=subset[,4]-subset[,5],
           subAminB=subset[,3]-subset[,4],
           subCminA=subset[,5]-subset[,3],
           subCminB=subset[,5]-subset[,4]) %>%
    select(1,2,starts_with("sub"))
  
  colnames(sbtrct)[3:7] = c(paste(genotype,tissues[1],"minus",tissues[3],sep="_"),
                            paste(genotype,tissues[2],"minus",tissues[3],sep="_"),
                            paste(genotype,tissues[1],"minus",tissues[2],sep="_"),
                            paste(genotype,tissues[3],"minus",tissues[1],sep="_"),
                            paste(genotype,tissues[3],"minus",tissues[2],sep="_"))
  
  
  # export the subtracted matrices
  for (i in 1:length(sbtrct)){
    a = i+2
    
    if (a > length(sbtrct)){
      break
    } else {
      
      export = 
        sbtrct %>%
        dplyr::select(1,2,a) %>%
        filter(.[,3]!=0)
      
      exp_file = paste0("Tiled-C_matNorm_subtract_", colnames(export)[3], ".matrix")
      
      write.table(export,
                  file=paste0(matrixFolder,exp_file),
                  quote=F,
                  row.names=F,
                  col.names=F)
    }
  }
}



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
for (i in 1:length(normed)){
  a = i+2
  
  if (a > length(normed)){
    break
  } else {
    
    export = 
      filtData %>%
      dplyr::select(starts_with("rel"),a) %>%
      filter(.[,3]!=0)
    
    exp_file = paste0("Tiled-C_matNorm_window_", chr, "-", UCSCmin,"-", UCSCmax,"_", binNum, "_bins_", colnames(export)[3], ".matrix")
    
    write.table(export,
                file=paste0(matrixFolder,exp_file),
                quote=F,
                row.names=F,
                col.names=F)
  }
}

  
cat("Export of window filtered normalised matrix for viewing window", minFrag, "to", maxFrag, "comlete", "\n")



###########################################
### FILTER GENO SUBTRACT DATA ON WINDOW ###
###########################################
# redo the GENOTYPE subtracts and then export the windowed matrices




for (tissue in tissues){
  
  # filter the union or normed table on just the desired minimum and maximum fragments
  filtData <- 
    normed %>%
    filter(baitID > minFrag &
             baitID < maxFrag &
             preyID > minFrag &
             preyID < maxFrag)
  
  
  subset =
    filtData %>%
    select(1,2,starts_with(tissue))
  
  sbtrct <-
    subset %>%
    mutate(subAminC=subset[,3]-subset[,5],
           subBminC=subset[,4]-subset[,5],
           subAminB=subset[,3]-subset[,4],
           subCminA=subset[,5]-subset[,3],
           subCminB=subset[,5]-subset[,4]) %>%
    select(1,2,starts_with("sub"))
  
  colnames(sbtrct)[3:7] = c(paste(tissue,genotypes[1],"minus",genotypes[3],sep="_"),
                            paste(tissue,genotypes[2],"minus",genotypes[3],sep="_"),
                            paste(tissue,genotypes[1],"minus",genotypes[2],sep="_"),
                            paste(tissue,genotypes[3],"minus",genotypes[1],sep="_"),
                            paste(tissue,genotypes[3],"minus",genotypes[2],sep="_"))
  
  
  
  # make into integers to allow filtering
  sbtrct$baitID <- as.integer(sbtrct$baitID)
  sbtrct$preyID <- as.integer(sbtrct$preyID)
  
  cat("Now selecting just a window of the data for tissue ", tissue, "\n")
  
  # for plotting in python, I need the bins numbered 0 onwards!
  # to get this just subtract the minFrag from each
  filtData <- 
    sbtrct %>%
    mutate(relBaitID=baitID-minFrag-1,
           relPreyID=preyID-minFrag-1)
  
  # rearrange filtData so relBaitID/ relPreyID are as cols 1 and 2 (needed for export of matrices)
  filtData %<>%
    select(starts_with("rel"),3:7)
  
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
  
  
  # export the subtracted matrices
  for (i in 1:length(filtData)){
    a = i+2
    
    if (a > length(filtData)){
      break
    } else {
      
      export = 
        filtData %>%
        dplyr::select(1,2,a) %>%
        filter(.[,3]!=0)
      
      exp_file = paste0("Tiled-C_matNorm_subtract_window_", chr, "-", UCSCmin,"-", UCSCmax,"_", binNum, "_bins_", colnames(export)[3], ".matrix")
      write.table(export,
                  file=paste0(matrixFolder,exp_file),
                  quote=F,
                  row.names=F,
                  col.names=F)
    }
  }
}


##############################################
### FILTER TISSUES SUBTRACT DATA ON WINDOW ###
##############################################
# redo the TISSUE subtracts and then export the windowed matrices


for (genotype in genotypes){
  
  # filter the union or normed table on just the desired minimum and maximum fragments
  filtData <- 
    normed %>%
    filter(baitID > minFrag &
             baitID < maxFrag &
             preyID > minFrag &
             preyID < maxFrag)
  
  
  subset =
    filtData %>%
    select(1,2,matches(genotype))
  
  sbtrct <-
    subset %>%
    mutate(subAminC=subset[,3]-subset[,5],
           subBminC=subset[,4]-subset[,5],
           subAminB=subset[,3]-subset[,4],
           subCminA=subset[,5]-subset[,3],
           subCminB=subset[,5]-subset[,4]) %>%
    select(1,2,starts_with("sub"))
  
  colnames(sbtrct)[3:7] = c(paste(genotype,tissues[1],"minus",tissues[3],sep="_"),
                            paste(genotype,tissues[2],"minus",tissues[3],sep="_"),
                            paste(genotype,tissues[1],"minus",tissues[2],sep="_"),
                            paste(genotype,tissues[3],"minus",tissues[1],sep="_"),
                            paste(genotype,tissues[3],"minus",tissues[2],sep="_"))
  
  
  # make into integers to allow filtering
  sbtrct$baitID <- as.integer(sbtrct$baitID)
  sbtrct$preyID <- as.integer(sbtrct$preyID)
  
  cat("Now selecting just a window of the data for tissue ", genotype, "\n")
  

  # for plotting in python, I need the bins numbered 0 onwards!
  # to get this just subtract the minFrag from each
  filtData <- 
    sbtrct %>%
    mutate(relBaitID=baitID-minFrag-1,
           relPreyID=preyID-minFrag-1)
  
  # rearrange filtData so relBaitID/ relPreyID are as cols 1 and 2 (needed for export of matrices)
  filtData %<>%
    select(starts_with("rel"),3:7)
  
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
  
  
  # export the subtracted matrices
  for (i in 1:length(filtData)){
    a = i+2
    
    if (a > length(filtData)){
      break
    } else {
      
      export = 
        filtData %>%
        dplyr::select(1,2,a) %>%
        filter(.[,3]!=0)
      
      exp_file = paste0("Tiled-C_matNorm_subtract_window_", chr, "-", UCSCmin,"-", UCSCmax,"_", binNum, "_bins_", colnames(export)[3], ".matrix")
      write.table(export,
                  file=paste0(matrixFolder,exp_file),
                  quote=F,
                  row.names=F,
                  col.names=F)
    }
  }
}



cat("Export of window filtered raw matrix for viewing window", minFrag, "to", maxFrag, "comlete", "\n")





#################################################
### FILTER NORMALISED DATA ON DESEQ DIFF INTS ###
#################################################
# as in a only keeping bins that have a significant interaction
# want to plot the normed reads

# load back the normalised reads calculated as above
#normed <- readRDS(paste0(variablesFolder, "/Tiled_mergeNormAllData.rds"))






# get the samples from the DESeq data directory
# only keep the text files with the actual comparisons!
comparisons <-
  list.files(path=outDirData) %>%
  data.frame %>%
  filter(grepl(".txt", .))

colnames(comparisons) <- "comparison"

comparisons = as.character(comparisons$comparison)

# just the first one
#comp=comparisons[1]


for (comp in comparisons){
  
  # load the DEseq2 table
  comp %<>% as.character
  compTable = read_tsv(paste0(outDirData,comp))
  
  if (nrow(compTable)==0){
    next
  }
  
  
  # now get the name of the samples in the comparison
  comp %<>%
    data.frame
  
  colnames(comp) <- "comp"
  
  comp %<>%
    separate(comp, into=c(NA,"cellType_A","genoType_A",NA,"cellType_B", "genoType_B"), sep="_") %>%
    mutate(sample_A=paste(cellType_A,genoType_A,"comb_normed",sep="_"),
           sample_B=paste(cellType_B,genoType_B,"comb_normed",sep="_"))
  
  
  sample_A=comp$sample_A
  sample_B=comp$sample_B
  
  # extract just these samples from the table
  
  normedA <- 
    normed %>%
    dplyr::select(1,2,matches(sample_A))
  
  normedB <- 
    normed %>%
    dplyr::select(1,2,matches(sample_B))
  
  # prep for join with DESeq frags
  toJoinA <-
    normedA %>%
    mutate(combo=paste0(baitID,"-",preyID)) %>%
    dplyr::select(combo,3)
  
  toJoinB <-
    normedB %>%
    mutate(combo=paste0(baitID,"-",preyID)) %>%
    dplyr::select(combo,3)
  
  
  # only keep the sig int frags and use log FC to chose which sample interaction is higher in
  diffIntA <- compTable %>% filter(padj<=FDRLevel&log2FoldChange>0)
  diffIntB <- compTable %>% filter(padj<=FDRLevel&log2FoldChange<0)
  
  
  diffIntA <-
    diffIntA %>%
    dplyr::select(frags,logpadj) %>%
    dplyr::rename(combo=frags)
  
  diffIntB <-
    diffIntB %>%
    dplyr::select(frags,logpadj) %>%
    dplyr::rename(combo=frags)
  
  
  
  
  # joining the normed data to the sigInt data and then only keeping the frags with values
  # splitting back up again
  # this should help plotting as I know exporting this df worked above
  joinedA <- 
    toJoinA %>%
    left_join(diffIntA,by="combo") %>%
    filter(logpadj>0) %>%
    separate(combo, into=c("baitID","preyID"), sep = "-") %>%
    dplyr::select(baitID,preyID,3)
  
  joinedB <- 
    toJoinB %>%
    left_join(diffIntB,by="combo") %>%
    filter(logpadj>0) %>%
    separate(combo, into=c("baitID","preyID"), sep = "-") %>%
    dplyr::select(baitID,preyID,3)
  
  
  
  
  
  ## subset the data on the desired window
  ## make into integers to allow filtering
  joinedA$baitID <- as.integer(joinedA$baitID)
  joinedA$preyID <- as.integer(joinedA$preyID)
  joinedB$baitID <- as.integer(joinedB$baitID)
  joinedB$preyID <- as.integer(joinedB$preyID)
  
  
  ### finding out the viewing window to use on UCSC
  # and get the binNum
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
  
  
  
  
  
  ## filter the normalised reads diff ints over the window
  
  joined_A_window <- 
    joinedA %>%
    filter(baitID > minFrag &
             baitID < maxFrag &
             preyID > minFrag &
             preyID < maxFrag)
  
  joined_B_window <- 
    joinedB %>%
    filter(baitID > minFrag &
             baitID < maxFrag &
             preyID > minFrag &
             preyID < maxFrag)
  
  # for plotting in python, I need the bins numbered 0 onwards!
  # to get this just subtract the minFrag from each
  # these are then relative to minFrag=0
  # -1 to correct for 0 indexing
  joined_A_window <- 
    joined_A_window %>%
    mutate(relBaitID=baitID-minFrag-1,
           relPreyID=preyID-minFrag-1) %>%
    dplyr::select(relBaitID,relPreyID,3)
  
  joined_B_window <- 
    joined_B_window %>%
    mutate(relBaitID=baitID-minFrag-1,
           relPreyID=preyID-minFrag-1) %>%
    dplyr::select(relBaitID,relPreyID,3)
  
  
  
  # make samples names back to P1-CTCF-KO (with - delimiter within genotype)
  # rather than P1.CTCF.KO (with . delimiter within genotype)
  # as all other exports have - delimiter
  # need to join sample names when plotting with fixed thresholds in python
  # so need to change them here and export as - delimiter
  
  sample_A = gsub("\\.", "-", sample_A)
  sample_B = gsub("\\.", "-", sample_B)
  
  
  ## only export the tables if they have any data!!
  if (nrow(joinedA)!=0){
    # export the DESeq filtered normalised reads significant matrices
    write.table(joinedA,
                file=paste0(matrixFolder,
                            "Tiled-C_matNorm_DESeq2_FDR_",
                            FDRLevel, "_",
                            paste(sample_A, "vs", sample_B,sep="_"),
                            "normReads.matrix"),
                quote=F,
                row.names=F,
                col.names=F)
  }
  
  if (nrow(joinedB)!=0){
    
    write.table(joinedB,
                file=paste0(matrixFolder,
                            "Tiled-C_matNorm_DESeq2_FDR_",
                            FDRLevel, "_",
                            paste(sample_B, "vs", sample_A,sep="_"),
                            "normReads.matrix"),
                quote=F,
                row.names=F,
                col.names=F)
  }
  

  
  
  
  if (nrow(joined_A_window)!=0){
    # export the DESeq filtered normalised reads and windowed significant matrices
    write.table(joined_A_window,
                file=paste0(matrixFolder,
                            "Tiled-C_matNorm_DESeq2_FDR_",
                            FDRLevel, "_",
                            "window_bins_",
                            binNum, "_", minFrag, "_to_", maxFrag, "_",
                            paste(sample_A, "vs", sample_B,sep="_"), "_",
                            "FDR_", FDRLevel, "normReads.matrix"),
                quote=F,
                row.names=F,
                col.names=F)
  }
  
  if (nrow(joined_B_window)!=0){
    
    write.table(joined_B_window,
                file=paste0(matrixFolder,
                            "Tiled-C_matNorm_DESeq2_FDR_",
                            FDRLevel, "_",
                            "window_bins_",
                            binNum, "_", minFrag, "_to_", maxFrag, "_",
                            paste(sample_B, "vs", sample_A,sep="_"), "_",
                            "FDR_", FDRLevel, "normReads.matrix"),
                quote=F,
                row.names=F,
                col.names=F)
  }
  
  cat("Export of window filtered normalised matrix for viewing window", minFrag, "to", maxFrag, "comlete for samples",sample_A, sample_B, "\n")
  
}



