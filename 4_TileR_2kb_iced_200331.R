### TILED ANALYSIS ###

library(tidyverse)
library(tools)

#### taking iced interaction matrices calculated by Marieke's script 
# from individual samples

# generating virtual Capture-C plots
# performing subtracts between tissues and between genotypes


# autohub_R at the end of this script will automatically generate the necessary hub files
# however it will append onto the text files if they exist
# so if an error occurs, make sure to delete the hub files before rerunning this script


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



###############
### OUTPUTS ###
###############

trackFolder <- paste0(base, "trackFolder_", matrixType, "/")
dir.create(trackFolder, showWarnings = F)

variablesFolder <- paste0(base, "savedVariables_", matrixType, "/")
dir.create(variablesFolder, showWarnings = F)

matrixFolder <- paste0(base, "outMatrixFolder_", matrixType, "/")
dir.create(matrixFolder, showWarnings = F)

## autohub.R
localhubDir <- paste0(base, "hubs_", matrixType, "/")
dir.create(localhubDir, showWarnings=F)



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
# exchange the - for . to prevent problems with column naming
genotypes <- gsub("-", ".", genotypes)
# clone names (reps effectively)
clones <- levels(factor(word(samples$sample, 3, sep="_")))





##################################
#### LOAD PREVIOUS STEP DATA #####
##################################

# iced individual replicate data was previously anlaysed in 3_TileR_replicateAnalysis_iced_normMatrices_200327
union <- readRDS(paste0(variablesFolder, "/Tiled_indiReps_", matrixType, ".rds"))

# make all columns of union numeric

for (i in 1:length(union)){
  union[,i] <- as.numeric(union[,i])
}



##########################
### BAIT SPECIFIC PART ###
##########################

# note:: this should be done on ICE normalised counts!!!! This script will then scaled to contact number as well.
# ICE normalisation does not scale for contact number between samples

# this part outputs virtual capC plots for each bait for each cell type / geno type group
# also performs subtracts between genotypes within each cell type
# and performs subtracts between cell types within each genotype



for (i in 1:nrow(baits)){ # working on each bait at a time
  
  # get the bait_frag, baitName, and baitChrom we are working on
  baitUsed = baits$bait_frag[i]
  baitName = baits$baitName[i]
  baitNum = baits$baitID[i]
  baitChrom = 
    baits[i,] %>%
    separate(bait_frag, into=c("baitChrom", NA), sep=":", remove=T) %>%
    dplyr::select(baitChrom) %>%
    as.character()
  
  cat("Now analysing bait", baits$baitName[i],  "Tile-C version", "\n")
  
  # filter the union table on just this bait
  # note I am actually not really using the concept of bait and prey
  # just looking for any fragment that interacts with this viewpoint!!
  # grep on the combo of "baitID-preyID" wasn't working as grep(187) found 1871 etc
  # also need to get a fragID column to plot the data over
  # also excluding the baited bin and +/- one either side
  
  countData <- 
    union %>%
    filter(baitNum==baitID | baitNum==preyID) %>%
    mutate(fragID=((baitID+preyID)-baitNum)) %>%
    filter(fragID!=baitNum) %>%
    filter(fragID!=baitNum+1) %>%
    filter(fragID!=baitNum-1) %>%
    dplyr::select(-baitID,-preyID) %>%
    tibble::column_to_rownames(var = "fragID")
  
  
  #####################
  ### NORMALISATION ###
  #####################
  
  ## exact normalisation commands taken from CC Stats pipeline ##
  ## put into the function above
  # hashed out the ones I don't want and added upper lower confidence
  #### this type of normalisation for Tiled data might be not appropriate, will try it and see
  ### if its no good then I will use raw counts matrix for DEseq, and ICE normalised matrix for plotting (this might be best anyway)
  
  # make into numeric to allow summing
  
  for (i in 1:length(countData)){
    countData[,i] <- as.numeric(countData[,i])
  }
  
  data.matrix <- data.matrix(countData) 
  column.totals <- apply(countData, 2, sum)
  data.norm <- t(t(data.matrix) * 5e3 / column.totals) 
  norm.totals <- apply(data.norm, 2, sum)
  column.totals
  norm.totals
  
  
  #cat("Performing the mean and sd for each tissue \n")
  # do the mean dif and sd for each group (combo of tissue and genotype in samples table)
  
  for (i in 1:nrow(samples)){
    
    group = paste(samples$tissue[i], samples$genotype[i],sep="_")
    
    # losing "-" in column names when passing data.norm to data.frame below
    group = gsub("-", ".", group)
    
    # get just the samples for this group from the normalised table
    group.data <-
      data.norm %>%
      data.frame %>%
      dplyr::select(matches(group))
      
    
    # get sum, mean, and sd (will have to specify diffs later as there could be many (by tissue and by genotype))
    group_mean <- data.frame(row.names=row.names(group.data),
                             apply(group.data, 1, mean))
    colnames(group_mean) <- paste0(group, "_mean")
                             
                             
    group_sd <- data.frame(row.names=row.names(group.data),
                           apply(group.data, 1, sd))
    colnames(group_sd) <- paste0(group, "_sd")
    
    #group_sum <- data.frame(row.names=row.names(group.data),
    #                        apply(group.data, 1, sum))
    
    # get upper confidence interval (mean+sd)
    group_upper <- group_mean + group_sd
    colnames(group_upper) <- paste0(group, "_upper_conf")
    
    
    # prepare a table for export
    #cat("Calculations done, now preparing files... \n")
    
    
    # make a list of all these dfs 
    # to allow formatting the bed coordinates
    # and exporting them in the below loop
    test <- list(group_mean,group_upper) # add others here
    
    #cat("Prep done, now exporting the files... \n")
    
    for (df in 1:length(test)){ # for each file export loop
      
      expo <- data.frame(test[df])
      fileName <- colnames(expo)
      
      # get the names of the fragments and convert to bedGraph style
      # convert 0 to missing data to save space in files
      expor <- expo %>% 
        mutate(fragID=as.integer(row.names(expo))) %>%
        left_join(fragData, by="fragID") %>%
        dplyr::select(prey_chr,prey_start,prey_end,1) %>%
        filter(.[,4]!=0)
      
      
      write.table(expor, file=paste0(trackFolder, baitName, "_", fileName, ".bdg"),
                  quote=F,
                  row.names = F,
                  col.names = F,
                  sep="\t")
      
    } # end of each file loop
    
    #cat("Export of normalised bedGraphs comlete for", baitName,  "\n")
    
  } # end of for sample group loop
  
  
  # now do the subtracts between genotypes within one tissue for this bait
  for (tissue in tissues){
    
    # get just the columns for this tissue from the table
    tissue.data <-
      data.norm %>%
      data.frame %>%
      dplyr::select(matches(tissue)) %>%
      rownames_to_column(var="intFrag")
    
    # change the genotype names to allow extracting from the table with - subbed for .
    genotypes2 = gsub("-", ".", genotypes)
  
    
    # convert to long format to make easier to group by and process means
    data_long <- gather(tissue.data, genotype, intValue, 2:ncol(tissue.data), factor_key=TRUE)
    

    # set up a genotype column
    data_long %<>%
      mutate(geno_group = if_else(condition = grepl(genotypes2[1], genotype), 
                              true=genotypes2[1],
                              if_else(condition = grepl(genotypes2[2], genotype),
                                      true=genotypes2[2],
                                      if_else(condition = grepl(genotypes2[3], genotype),
                                                        true=genotypes2[3], false="none"))))
      
    
    
    # get the mean interactions in each intFrag by genotype
    data_long %<>%
    group_by(geno_group,intFrag) %>%
      mutate(group_mean=mean(intValue))
    
    # make a new df which only keeps one intFrag for each group 
    # (loses replicate info)
    combined <-
    data_long %>%
      dplyr::select(intFrag,geno_group,group_mean) %>%
      filter(row_number(intFrag) == 1)
    
    
    P1 <-
      combined %>%
      filter(geno_group==genotypes2[1])
    
    P2 <-
      combined %>%
      filter(geno_group==genotypes2[2])
    
    WT <-
      combined %>%
      filter(geno_group==genotypes2[3])
    
    
    # combine into wide format and then do subtractions
    tmp <-
      P1 %>% ungroup %>%
      left_join(P2, by="intFrag") %>%
      left_join(WT, by="intFrag") %>%
      dplyr::select(intFrag,matches("mean")) %>%
      dplyr::rename(P1.CTCF.KO_mean=group_mean.x,
                    P2.CTCF.KO_mean=group_mean.y,
                    WT_mean=group_mean) %>%
      mutate(WT_minus_P1.CTCF.KO=WT_mean-P1.CTCF.KO_mean,
             WT_minus_P2.CTCF.KO=WT_mean-P2.CTCF.KO_mean,
             P1.CTCF.KO_minus_P2.CTCF.KO=P1.CTCF.KO_mean-P2.CTCF.KO_mean) %>%
      mutate(intFrag=as.numeric(intFrag)) %>%
      arrange(intFrag) %>%
      tibble::column_to_rownames(var = "intFrag") %>%
      dplyr::select(matches("minus"))
    
    
    # now export each of the columns of tmp 
    # in the loop below
    # and format the bed coordinates

    
    for (i in 1:ncol(tmp)){ # for each file export loop
      
      expo <- data.frame(row.names=row.names(tmp), data=tmp[,i])
      colnames(expo) <- paste0(tissue, "_", colnames(tmp)[i])
      fileName <- colnames(expo)

      # get the names of the fragments and convert to bedGraph style
      # convert 0 to missing data to save space in files
      expor <- expo %>% 
        mutate(fragID=as.integer(row.names(expo))) %>%
        left_join(fragData, by="fragID") %>%
        dplyr::select(prey_chr,prey_start,prey_end,1) %>%
        filter(.[,4]!=0)
      
      
      write.table(expor, file=paste0(trackFolder, baitName, "_", fileName, ".bdg"),
                  quote=F,
                  row.names = F,
                  col.names = F,
                  sep="\t")
      
    } # end of each file loop
    
    
  } # end of for tissues loop (within one bait)
  

  
  # now do the subtracts between tissues within one genotype for this bait
  for (genotype in genotypes2){
    
    # get just the columns for this tissue from the table
    tissue.data <-
      data.norm %>%
      data.frame %>%
      dplyr::select(matches(genotype)) %>%
      rownames_to_column(var="intFrag")
    
    
    # convert to long format to make easier to group by and process means
    data_long <- gather(tissue.data, genotype, intValue, 2:ncol(tissue.data), factor_key=TRUE)
    
    
    # set up a tissue column
    data_long %<>%
      mutate(tissue_group = if_else(condition = grepl(tissues[1], genotype), 
                                  true=tissues[1],
                                  if_else(condition = grepl(tissues[2], genotype),
                                          true=tissues[2],
                                          if_else(condition = grepl(tissues[3], genotype),
                                                  true=tissues[3], false="none"))))
    
    
    
    # get the mean interactions in each intFrag by genotype
    data_long %<>%
      group_by(tissue_group,intFrag) %>%
      mutate(group_mean=mean(intValue))
    
    # make a new df which only keeps one intFrag for each group 
    # (loses replicate info)
    combined <-
      data_long %>%
      dplyr::select(intFrag,tissue_group,group_mean) %>%
      filter(row_number(intFrag) == 1)
    
    
    CD41 <-
      combined %>%
      filter(tissue_group==tissues[1])
    
    Flk1 <-
      combined %>%
      filter(tissue_group==tissues[2])
    
    Un <-
      combined %>%
      filter(tissue_group==tissues[3])
    

    # combine into wide format
    # then do subtractions
    # then only keep the subtraction columns as means and upper confs have already been calculated
    tmp2 <-
      CD41 %>% ungroup %>%
      left_join(Flk1, by="intFrag") %>%
      left_join(Un, by="intFrag") %>%
      dplyr::select(intFrag,matches("mean")) %>%
      dplyr::rename(CD41_mean=group_mean.x,
                    Flk1_mean=group_mean.y,
                    Un_mean=group_mean) %>%
      mutate(CD41_minus_Flk1=CD41_mean-Flk1_mean,
             CD41_minus_Un=CD41_mean-Un_mean,
             Flk1_minus_Un=Flk1_mean-Un_mean) %>%
      mutate(intFrag=as.numeric(intFrag)) %>%
      arrange(intFrag) %>%
      tibble::column_to_rownames(var = "intFrag") %>%
      dplyr::select(matches("minus"))
      

    
    for (i in 1:ncol(tmp2)){ # for each file export loop
      
      expo <- data.frame(row.names=row.names(tmp2), data=tmp2[,i])
      colnames(expo) <- paste0(genotype, "_", colnames(tmp2)[i])
      fileName <- colnames(expo)
      
      # get the names of the fragments and convert to bedGraph style
      # convert 0 to missing data to save space in files
      expor <- expo %>% 
        mutate(fragID=as.integer(row.names(expo))) %>%
        left_join(fragData, by="fragID") %>%
        dplyr::select(prey_chr,prey_start,prey_end,1) %>%
        filter(.[,4]!=0)
      
      
      write.table(expor, file=paste0(trackFolder, baitName, "_", fileName, ".bdg"),
                  quote=F,
                  row.names = F,
                  col.names = F,
                  sep="\t")
      
    } # end of each file loop
    
  
  } # end of for each genotype loop
  
  

  
} # end of for bait loop



####################
##### autoHub.R ####
####################


# this script takes the list of baits, and the tissue types from the data folder, 
# and outputs the details of UCSC hub that should be placed into 
# publicDir eg (/path/to/public/hub/asSpecified/) by the user

# each bait for each cell type gets these tracks :

# overlay track of mean and mean + standard deviation (upper_conf)

#### make sure to copy over the files from:

# base/trackFolder

#
# then convert to bigWig
# module load ucsctools
# fetchChromSizes mm9 > mm9.chrom.sizes
# for file in * ; do bedSort $file $file ; bedGraphToBigWig $file mm9.chrom.sizes "${file}.bw" ; done

# Then move to:
# publicDir

cat("Now building the UCSC hub... \n")


# first set up a hub one for each group (cellType_genoType)

samples2 <-
  samples %>%
    mutate(group=paste0(tissue,"_",genotype)) %>%
    dplyr::select(group) %>% unique

groups = samples2$group
# losing "-" in column names when passing data.norm to data.frame below
groups = gsub("-", ".", groups)


for (a in 1:length(groups)) { ### begin group loop
  
  # get the cell type
  # create new files in a cellType specific directory 
  # one directory / hub for each cell type
  cellType <- groups[a]
  tissuehubDir <- paste0(localhubDir, "hub_", cellType, "_", matrixType, "/")
  dir.create(tissuehubDir, showWarnings = F)
  hubDir <- paste0(publicDir, "hub_", cellType, "_", matrixType, "/") 
  
  ## match colours to the different groups
  
  # define each cell types track colours a= green themed, b=blue themed like Jelena's CC / stats pipes
  if (grepl("CD41",cellType)) {
    mainCol <- "112,179,143" # darkgreen
    subCol <- "205,222,216" # lightgreen
    #subCol <- "255,255,255" # white
  } else if (grepl("Un",cellType)) {
    mainCol <- "133,149,194" # darkblue
    subCol <- "183,195,208" # lightblue
    #subColLow <- "255,255,255" # white
  } else if (grepl("Flk1",cellType)) {
    mainCol <- "201,178,211" # darkpurple
    subCol <- "252,231,255" # lightpurple
  }
  
  ## begin making the cell type specific hub ##
  localhubDirCell <- paste0(localhubDir, "hub_", cellType,"_", matrixType, "/")
  
  # tracks file written to in the oligos loop
  tracksFile <- paste0(localhubDirCell, cellType, "_", matrixType, "_tracks.txt")
  
  # hub file created now
  hubFile <- paste0(localhubDirCell, cellType, "_", matrixType, "_hub.txt")
  hubName <- paste0("autoHubR_", cellType, "_", matrixType, "2kbRes_Virtual_CapC")
  
  cat("hub ", hubName, "\n",
      "shortLabel ", hubName, "\n",
      "longLabel ", hubName, "\n",
      "genomesFile ", hubDir, cellType, "_", matrixType, "_genomes.txt", "\n",
      "email ", email,
      file=hubFile,
      append=F,
      sep="")
  
  # genomes file
  genomesFile <- paste0(localhubDirCell, cellType, "_", matrixType, "_genomes.txt")
  cat("genome ", genome, "\n",
      "trackDb ", hubDir, cellType, "_", matrixType, "_tracks.txt", "\n",
      file=genomesFile,
      append=F,
      sep="")
  
  # hubAddress file
  hubAddressFile <- paste0(localhubDirCell, cellType, "_", matrixType, "_hubAddress.txt")
  cat(hubDir, cellType, "_", matrixType, "_hub.txt",
      file=hubAddressFile,
      append=F,
      sep="")
  
  ## loop through each bait and append to tracks file
  
  for (r in 1:nrow(baits)) {  ### begin oligo loop
    
    thisBait <- baits[[r,3]]
    
    # variables updated within the oligo loop
    ## the names of all the tracks to include per bait
    
    ## overlay track
    trackName <- paste(thisBait, cellType, "virtCapC_meanOverlay", sep="_")
    meanName <- paste(thisBait, cellType, "mean", sep="_")
    upperName <- paste(thisBait, cellType, "upper_conf", sep="_")

    
    header <- paste(paste0(rep("#",10),collapse=""), cellType, matrixType, thisBait, paste0(rep("#",20),collapse=""), sep="   ")
    
    # overlay track
    cat("\n", "\n", header, "\n", "\n",  # main tracks.txt file cat command start
        "track ", trackName, "\n",
        "type bigWig \n",
        "container multiWig", "\n",
        "shortLabel ", trackName, "\n",
        "longLabel ", trackName, "\n",
        "visibility full", "\n",
        "aggregate transparentOverlay", "\n",
        "showSubtrackColorOnUi on", "\n",
        "maxHeightPixels 100:50:0", "\n",
        #"autoScale on", "\n",
        "viewLimits 0:63", "\n",
        "windowingFunction mean", "\n",
        #"priority 1", "\n",
        "\n",  "\t",
        
        "track ", upperName, "\n", "\t",
        "type bigWig \n", "\t",
        "parent ", trackName, "\n", "\t",
        "color ", subCol, "\n", "\t",
        "bigDataUrl ", publicDir, upperName, ".bdg.bw", "\n", "\t",
        "shortLabel ", upperName, "\n", "\t",
        "longLabel ", upperName, "\n", "\t",
        #"autoScale on", "\n",
        "\n",  "\t",
        
        "track ", meanName, "\n", "\t",
        "type bigWig \n", "\t",
        "parent ", trackName, "\n", "\t",
        "color ", mainCol, "\n", "\t",
        "bigDataUrl ", publicDir, meanName, ".bdg.bw", "\n", "\t",
        "shortLabel ", meanName, "\n", "\t",
        "longLabel ", meanName, "\n", "\t",
        #"autoScale on", "\n",
        "\n",   "\t",
        

        file = tracksFile, # tracks file
        append=T, # append or overwrite
        sep="") # end of main cat command
    
    
  } ### end bait loop
  
  
  # now list the files in the track folder
  # get just a list of all the subtracts
  # for this group
  # add these to the tracks file for this group
  
  subs <-
    list.files(path=trackFolder) %>%
    data.frame %>%
    filter(grepl("minus", `.`) &
             !grepl("mean", `.`) &
             !grepl("upper_conf", `.`) &
             grepl(cellType, `.`)) %>%
    mutate(subtracts=file_path_sans_ext(`.`)) %>%
    dplyr::select(-`.`)
  
  subtracts=subs$subtracts
  
  
  header2 <- paste(paste0(rep("#",10),collapse=""), cellType, matrixType, "subtracts", paste0(rep("#",20),collapse=""), sep="   ")
  
  cat("\n", "\n", header2, "\n", "\n",
      file = tracksFile, # tracks file
      append=T, # append or overwrite
      sep="")
  
  
  # now do a for loop over each subtract to add it to the tracks file
  for (sub in subtracts){ # start of subtract file loop
    
    ## subtract track
    subTrackName <- sub
    subTrackFile <- paste0(sub, ".bdg.bw")

    # overlay track
    cat("track ", subTrackName, "\n",
      "type bigWig \n",
      "color 153,153,153", "\n",
      "maxHeightPixels 100:40:0", "\n",
      "visibility hide", "\n",
      "bigDataUrl ", publicDir, subTrackFile, "\n",
      "shortLabel ", subTrackName, "\n",
      "longLabel ", subTrackName, "\n",
      "autoScale on", "\n",
      "windowingFunction mean", "\n",
      "\n",  "\n",
        
      file = tracksFile, # tracks file
      append=T, # append or overwrite
      sep="") # end of main cat command
    
    
    
  } # end of subtract file loop

  
} ### end group loop (cellType_genoType)




# next set up hubs for each genotype
# losing "-" in column names when passing data.norm to data.frame below
samples$genotype = gsub("-", ".", samples$genotype)


samples2 <-
  samples %>%
  mutate(group=paste0(tissue,"_",genotype)) %>%
  dplyr::select(group) %>% unique

groups = samples2$group



for (geno in genotypes) { ### begin genotype loop
  
  # subset the samples table on just this genotype
  # then combine the tissue and genotype info to get group
  samples3 <-
    samples %>%
    filter(genotype==geno) %>%
    mutate(group=paste0(tissue, "_", genotype))
  
  
  genoGroups <- levels(factor(samples3$group))
  
  
  for (genoGroup in genoGroups){
    
    # get the cell type
    # create new files in a cellType specific directory 
    # one directory / hub for each cell type
    cellType <- geno
    tissuehubDir <- paste0(localhubDir, "hub_", cellType, "_", matrixType, "/")
    dir.create(tissuehubDir, showWarnings = F)
    hubDir <- paste0(publicDir, "hub_", cellType, "_", matrixType, "/") 
    
    ## match colours to the different groups
    
    # define each cell types track colours a= green themed, b=blue themed like Jelena's CC / stats pipes
    if (grepl("CD41",genoGroup)) {
      mainCol <- "112,179,143" # darkgreen
      subCol <- "205,222,216" # lightgreen
      #subCol <- "255,255,255" # white
    } else if (grepl("Un",genoGroup)) {
      mainCol <- "133,149,194" # darkblue
      subCol <- "183,195,208" # lightblue
      #subColLow <- "255,255,255" # white
    } else if (grepl("Flk1",genoGroup)) {
      mainCol <- "201,178,211" # darkpurple
      subCol <- "252,231,255" # lightpurple
    }
    
    ## begin making the genotype specific hub ##
    localhubDirCell <- paste0(localhubDir, "hub_", cellType, "_", matrixType, "/")
    
    # tracks file written to in the oligos loop
    tracksFile <- paste0(localhubDirCell, cellType, "_", matrixType, "_tracks.txt")
    
    # hub file created now
    hubFile <- paste0(localhubDirCell, cellType,  "_", matrixType, "_hub.txt")
    hubName <- paste0("autoHubR_", cellType, "_", matrixType, "2kbRes_Virtual_CapC")
    
    cat("hub ", hubName, "\n",
        "shortLabel ", hubName, "\n",
        "longLabel ", hubName, "\n",
        "genomesFile ", hubDir, cellType,  "_", matrixType, "_genomes.txt", "\n",
        "email ", email,
        file=hubFile,
        append=F,
        sep="")
    
    # genomes file
    genomesFile <- paste0(localhubDirCell, cellType,  "_", matrixType, "_genomes.txt")
    cat("genome ", genome, "\n",
        "trackDb ", hubDir, cellType,  "_", matrixType, "_tracks.txt", "\n",
        file=genomesFile,
        append=F,
        sep="")
    
    # hubAddress file
    hubAddressFile <- paste0(localhubDirCell, cellType,  "_", matrixType, "_hubAddress.txt")
    cat(hubDir, cellType,  "_", matrixType, "_hub.txt",
        file=hubAddressFile,
        append=F,
        sep="")
    
    ## loop through each bait and append to tracks file
    
    for (r in 1:nrow(baits)) {  ### begin oligo loop
      
      thisBait <- baits[[r,3]]
      
      # variables updated within the oligo loop
      ## the names of all the tracks to include per bait
      
      ## overlay track information
      # need to put a and b infront of tracks so that they are plotted in the correct order!
      trackName <- paste(thisBait, genoGroup, "virtCapC_meanOverlay", sep="_")
      meanName <- paste(thisBait, genoGroup, "mean", sep="_")
      upperName <- paste(thisBait, genoGroup, "upper_conf", sep="_")
      
      
      header <- paste(paste0(rep("#",10),collapse=""), genoGroup, matrixType, thisBait, paste0(rep("#",20),collapse=""), sep="   ")
       
      # overlay track
      cat("\n", "\n", header, "\n", "\n",  # main tracks.txt file cat command start
          "track ", trackName, "\n",
          "type bigWig \n",
          "container multiWig", "\n",
          "shortLabel ", trackName, "\n",
          "longLabel ", trackName, "\n",
          "visibility full", "\n",
          "aggregate transparentOverlay", "\n",
          "showSubtrackColorOnUi on", "\n",
          "maxHeightPixels 100:50:0", "\n",
          #"autoScale on", "\n",
          "viewLimits 0:63", "\n",
          "windowingFunction mean", "\n",
          #"priority 1", "\n",
          "\n",  "\t",
          
          "track ", upperName, "\n", "\t",
          "type bigWig \n", "\t",
          "parent ", trackName, "\n", "\t",
          "color ", subCol, "\n", "\t",
          "bigDataUrl ", publicDir, upperName, ".bdg.bw", "\n", "\t",
          "shortLabel ", upperName, "\n", "\t",
          "longLabel ", upperName, "\n", "\t",
          #"autoScale on", "\n",
          "\n",  "\t",
          
          "track ", meanName, "\n", "\t",
          "type bigWig \n", "\t",
          "parent ", trackName, "\n", "\t",
          "color ", mainCol, "\n", "\t",
          "bigDataUrl ", publicDir, meanName, ".bdg.bw", "\n", "\t",
          "shortLabel ", meanName, "\n", "\t",
          "longLabel ", meanName, "\n", "\t",
          #"autoScale on", "\n",
          "\n",   "\t",
          
          
          file = tracksFile, # tracks file
          append=T, # append or overwrite
          sep="") # end of main cat command
      
      
    } ### end bait loop
    
  } # end genoGroup (cellType_genotype)
  
  
  # now list the files in the track folder
  # get just a list of all the subtracts
  # for this group
  # add these to the tracks file for this group
  
  subs <-
    list.files(path=trackFolder) %>%
    data.frame %>%
    filter(grepl("minus", `.`) &
             !grepl("mean", `.`) &
             !grepl("upper_conf", `.`) &
             grepl(geno, `.`) &
             !grepl("WT_minus_P1",`.`) &
             !grepl("WT_minus_P2", `.`)) %>%
    mutate(subtracts=file_path_sans_ext(`.`)) %>%
    dplyr::select(-`.`)
  
  subtracts=subs$subtracts
  
  
  header2 <- paste(paste0(rep("#",10),collapse=""), cellType, matrixType, "subtracts", paste0(rep("#",20),collapse=""), sep="   ")
  
  cat("\n", "\n", header2, "\n", "\n",
      file = tracksFile, # tracks file
      append=T, # append or overwrite
      sep="")
  
  
  # now do a for loop over each subtract to add it to the tracks file
  for (sub in subtracts){ # start of subtract file loop
    
    ## subtract track
    subTrackName <- sub
    subTrackFile <- paste0(sub, ".bdg.bw")
    
    # overlay track
    cat("track ", subTrackName, "\n",
        "type bigWig \n",
        "color 153,153,153", "\n",
        "maxHeightPixels 100:40:0", "\n",
        "visibility hide", "\n",
        "bigDataUrl ", publicDir, subTrackFile, "\n",
        "shortLabel ", subTrackName, "\n",
        "longLabel ", subTrackName, "\n",
        "autoScale on", "\n",
        "windowingFunction mean", "\n",
        "\n",  "\n",
        
        file = tracksFile, # tracks file
        append=T, # append or overwrite
        sep="") # end of main cat command
    
    
    
  } # end of subtract file loop
  
  
  
} ### end genotype loop



