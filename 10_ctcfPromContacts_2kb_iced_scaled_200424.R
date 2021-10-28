### TILED ANALYSIS ###

library(tidyverse)
library(tools)
library(cowplot)
library(magrittr)
library(RColorBrewer)
library(broom)
library(plotrix)


#### taking iced and scaled interaction matrices calculated by Marieke's script 
# from individual samples

# calculating CTCF promoter contacts in different groups (with replicate info as SD)


# working xx xx xxxx


###############
### OPTIONS ###
###############

options(max.print=100)


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


## ctcfData
# had bed file of ctcf peak cordinates and using bedtools 
# intersected with fragData to get overlapping bins
ctcfFile <- paste0(base, "ctcfData_bins_clip.txt")


###############
### OUTPUTS ###
###############

outFolder <- paste0(base, "ctcfProm_contacts_", matrixType, "/")
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

# prepare the enhancer data table to allow for some enhancers spread over multiple bins

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


########################################
#### GATHER CTCF PROMOTER CONTACTS #####
########################################
totals = data.frame()

for (i in 1:nrow(promData)){ # for prom loop
  
  promName=promData$promName[i]
  promID=promData$promID[i]
  
  
  # filter the matrix on just this prom
  prom_filt <-
    d %>% 
    filter(baitID==promID | preyID==promID) %>%
    mutate(promName=promName)
  
  # filter out all the excluded delIDs
  for (delID in delIDs){
    prom_filt %<>% 
      filter(!baitID %in% delID &
               !preyID %in% delID)
  }

  
  # start the for enhancer loop

  for (e in 1:nrow(ctcfData_new)){
    
    # get the info for just this enhancer
    ctcf_numName=ctcfData_new$ctcf_numName[e]
    #ctcfName=ctcfData_new$ctcfName[e]
    ctcfID=ctcfData_new$ctcfID[e]
    
    
    # filter on just this ctcfancers interactions with just this promoter
    ctcf_filt <-
      prom_filt %>% 
      filter(baitID==ctcfID | 
               preyID==ctcfID) %>%
      mutate(ctcf_numName=ctcf_numName)
    
    # add to output df
    totals = rbind.data.frame(totals,ctcf_filt)
    
  } # end of for ctcfancer loop
  
} # end of for prom loop




# get total E-P contacts in each group
# start by getting it for each sample
# and then split up the sample into group
sum_CTCF_contacts <-
  totals %>%
  group_by(sample,promName) %>%
  summarise(sum=sum(contacts)) %>%
  separate(sample,into=c("tissue","genotype", "clone", "XP", "ID"),sep="_") %>%
  mutate(group=paste0(tissue,"_",genotype)) %>%
  ungroup()




###########################################
#### ANOVA ON TOTAL PROMOTER CONTACTS #####
###########################################
## okay results so far, working code

# could be best to extract the enhancers and promoters from the DESeq table

# do stats (Four way anova)
res.aov <- aov(contacts ~ tissue*genotype*promName*ctcf_numName, data = totals)
fourWay = tidy(res.aov)
# get the significant effects
sigEffects = fourWay[fourWay$p.value<0.05,]$term
sigEffects <- sigEffects[!is.na(sigEffects)]

# do the post-hocTest and keep only if that main / interaction effect is significant
fourWay_phoc=
  tidy(TukeyHSD(res.aov)) %>%
  filter(term %in% sigEffects) %>%
  filter(adj.p.value<0.1)


# export the ANOVA table
write.table(fourWay,
            file=paste0(outFolder,"Four-Way-ANOVA_Summary-tables.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            na = "",
            sep="\t")

# export the post Hoc results
write.table(fourWay_phoc,
            file=paste0(outFolder,"Four-Way-ANOVA_postHocTestingOnlySigEffects.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            na = "",
            sep="\t")

# Now doing smaller Two-Way ANOVAs on each genotype/promoter separately

out_tables <- data.frame()
out_postHoc <- data.frame()
out_postHoc_main <- data.frame()

for (geno in genotypes){
  
  for (i in 1:nrow(promData)){
    
    prom=promData$promName[i]
    
    # filter the table on just this geno
    totals_temp <-
      filter(totals, genotype==geno & promName==prom)
    
    # perform the ANOVA
    res.aov <- aov(contacts ~ tissue*ctcf_numName, data = totals_temp)
    
    # get the ANOVA table and add to main table
    # using tidy from broom package
    out <- broom::tidy(res.aov)
    # add header and tail
    out = head_tail(out,paste("Two-way-ANOVA-table",geno,prom,sep="_"))
    # get the terms as factor so the characters can be changed to numeric in next step
    out$term=factor(out$term)
    #convert values to numeric
    out %<>%
      mutate_if(is.character,as.numeric)
    # add to main df
    out_tables <- rbind.data.frame(out_tables, out)
    
    
    
    # do the post-hoc tests
    # posthoc test just to keep the main effect of tissue
    postHoc_tissue <-
      tidy(TukeyHSD(res.aov)) %>%
      filter(term=="tissue")
    # get the table ready to rbind with the others
    postHoc_tissue = head_tail(postHoc_tissue,paste("Two-way-ANOVA-tissue-mainEffect",geno,prom,sep="_"))
    # add to main df
    out_postHoc_main <- rbind.data.frame(out_postHoc_main, postHoc_tissue)
    
    
    
    # now post hoc tests to get the significant interaction effects
    # to find enhancers that interact more or less in different tissues
    postHoc_int <-
      tidy(TukeyHSD(res.aov)) %>%
      filter(term=="tissue:ctcf_numName" &
               adj.p.value<0.2) %>%
      separate(comparison, into=c("comp1","comp2"),sep="-") %>%
      separate(comp1, into=c("tissue1","ctcf1"),sep=":") %>%
      separate(comp2, into=c("tissue2","ctcf2"),sep=":") %>%
      filter(ctcf1==ctcf2)
  
    # get the table ready to rbind with the others
    postHoc_int = head_tail(postHoc_int,paste("Two-way-ANOVA-tissue-enh-intEffect",geno,prom,sep="_"))
    # get the terms as factor so the characters can be changed to numeric in next step
    postHoc_int$term=factor(postHoc_int$term)
    postHoc_int$tissue1=factor(postHoc_int$tissue1)
    postHoc_int$tissue2=factor(postHoc_int$tissue2)
    #convert values to numeric
    postHoc_int %<>%
      mutate_if(is.character,as.numeric)
    # add to main df
    out_postHoc <- rbind.data.frame(out_postHoc, postHoc_int)
    
  } # end prom loop
  
} # end geno loop


# export the ANOVA tables (one anova for each genotype)
write.table(out_tables,
            file=paste0(outFolder,"Two-Way-ANOVA-by-genotype-by-promoter_Summary-tables.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            na = "",
            sep="\t")

# export the posthoc results
write.table(out_postHoc,
            file=paste0(outFolder,"Two-Way-ANOVA-by-genotype-by-promoter_postHoc-results_intEffect.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            na = "",
            sep="\t")

write.table(out_postHoc_main,
            file=paste0(outFolder,"Two-Way-ANOVA-by-genotype-by-promoter_postHoc-results_mainEffect_Tissue.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            na = "",
            sep="\t")



###############################################
#### PLOT SEPARATE CTCF PROMOTER CONTACTS #####
###############################################


### NOT WORKING FOR CTCF FOR SOME REASON

# now calculate mean and sd of the sum of contacts for each group
# give the genotypes a shorter name
output <- 
  totals %>%
  group_by(tissue,genotype,promName,ctcf_numName) %>%
  summarise(mean=mean(contacts),
            sd=sd(contacts)) %>%
  ungroup() %>%
  mutate(genotype=str_replace(genotype,"\\.CTCF\\.KO", "KO"),
         promName=str_replace(promName,"Runx1\\.", ""))


# get the order of the cell types proper
output$tissue = factor(output$tissue, levels=c("Un", "Flk1", "CD41"))


# plot each enhancer separately
plt_list <- list()
for (i in 1:nrow(ctcfData_new)){ # for enh_numName loop
  
  # get the enhancer names
  ctcf=ctcfData_new$ctcf_numName[i]
  
  temp <-
  output %>%
    filter(ctcf_numName==ctcf)
  
  
  # plot these by cell type, and promoter
  plt_list[[ctcf]] <-
    ggplot(temp, aes(x=tissue,y=mean))+
    geom_errorbar(aes(ymin=mean+sd, ymax=mean+sd),position=position_dodge(.9),width=.2)+
    geom_linerange(aes(ymin = mean, ymax = mean+sd))+
    geom_bar(stat = "identity",position=position_dodge(), aes(fill=tissue))+
    scale_fill_brewer(palette="Accent")+
    facet_grid(promName~genotype+ctcf_numName)+
    ggtitle(ctcf)+
    guides(fill=FALSE)+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid = element_blank(),
          strip.background =element_rect(fill="white"))
  

  
} # end of for ctcf_numName loop

  

plt_all <- plot_grid(plotlist=plt_list,ncol=6)


ggsave(plot=plt_all,
       filename=paste0(outFolder, "ctcfPromContacts_all.pdf"),
       width=12,
       height = 12)


#################################
#### PLOT A SUBSET OF CTCFs #####
#################################

# just using this to plot all CTCFs!!!
good_ctcf = levels(factor(output$ctcf_numName))

# or actually plotting a subset
#good_ctcf = c(-628.8,-467.6,-311.8,259.8,281.2,305.6)


# filter the already calcualted mean contacts and sd in different cell type / genotype groups
ctcf_filt <-
  filter(output, ctcf_numName %in% good_ctcf)

# filter the enhancer data table as well
ctcfData_filt <-
  filter(ctcfData_new, ctcf_numName %in% good_ctcf)


# get the order of the cell types proper
ctcf_filt$tissue = factor(ctcf_filt$tissue, levels=c("Un", "Flk1", "CD41"))

# plot each ctcfancer separately
plt_list <- list()
for (i in 1:nrow(ctcfData_filt)){ # for ctcf_numName loop
  
  # get the ctcfancer names
  ctcf=ctcfData_filt$ctcf_numName[i]
  
  temp <-
    ctcf_filt %>%
    filter(ctcf_numName==ctcf)
  
  for (i in 1:nrow(promData)){ # for prom loop
    
    # get the promoter names
    prom=strsplit(promData$promName[i],"\\.")[[1]][2]
    promID=promData$promID[i]
    
    ctcf_prom_filt <-
      temp %>%
      filter(promName==prom)
    
    # get a name for this plot combo
    combo = paste0(ctcf," (",prom,")")
    
    # plot these by cell type, and promoter
    plt_list[[combo]] <-
      ggplot(ctcf_prom_filt, aes(x=tissue,y=mean))+
      geom_errorbar(aes(ymin=mean+sd, ymax=mean+sd),position=position_dodge(.9),width=.2)+
      geom_linerange(aes(ymin = mean, ymax = mean+sd))+
      geom_bar(stat = "identity",position=position_dodge(), aes(fill=tissue))+
      scale_fill_brewer(palette="Accent")+
      facet_grid(promName~genotype)+
      ggtitle(combo)+
      guides(fill=FALSE)+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            panel.grid = element_blank(),
            strip.background =element_rect(fill="white"))
    
    
  } # end of for prom loop
  
  
} # end of for ctcf_numName loop


plt_all <- plot_grid(plotlist=plt_list,ncol=6)


ggsave(plot=plt_all,
       filename=paste0(outFolder, "ctcfPromContacts_all.pdf"),
       width=20,
       height = 20)

############################################
#### PLOT TOTAL CTCF PROMOTER CONTACTS #####
############################################

total_EPs <-
  totals %>%
  group_by(tissue,genotype,promName) %>%
  summarise(mean=mean(contacts),
            sd=std.error(contacts))

# get the order of the cell types proper
total_EPs$tissue = factor(total_EPs$tissue, levels=c("Un", "Flk1", "CD41"))


plot_EPs <-
  ggplot(total_EPs, aes(x=tissue,y=mean))+
  geom_errorbar(aes(ymin=mean+sd, ymax=mean+sd),position=position_dodge(.9),width=.2)+
  geom_linerange(aes(ymin = mean, ymax = mean+sd))+
  geom_bar(stat = "identity",position=position_dodge(), aes(fill=tissue))+
  scale_fill_brewer(palette="Accent")+
  facet_grid(promName~genotype)+
  ggtitle("Total CTCF-Promoter Contacts")+
  #guides(fill=T)+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        strip.background =element_rect(fill="white"))


ggsave(plot=plot_EPs,
       filename=paste0(outFolder, "ctcfPromContacts_total_se.pdf"),
       width=6,
       height = 3)