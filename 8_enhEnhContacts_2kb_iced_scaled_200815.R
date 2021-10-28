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

# calculating enhancer ENHANCER contacts in different groups (with replicate info as SD)

## will also want to subset on just the first intron of Runx1 as P1-CTCF-KO appears to alter this


# working xx xx 2020


## NOT FINISHED BUT NOTHING CLEAR IN RESULTS SO LEAVING INCOMPLETE


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


## enhData
# had bed file of enhancer cordinates and using bedtools 
# intersected with fragData to get overlapping bins
enhFile <- paste0(base, "enhData_bins.txt")


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



### just keep the enhancers in intron 1
enhData_new %<>%
  filter(start>92696373&stop<92826425 &
           enh_numName %in% c(110.2,24.2,3.2))


## now just keep the main enhancers we care about


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
#### GATHER ENHANCER ENHANCER CONTACTS #####
############################################
totals = data.frame()

for (i in 1:nrow(enhData_new)){ # for enh loop
  
  enhName=enhData_new$enhName[i]
  enhID=enhData_new$enhID[i]
  
  
  # filter the matrix on just this enh
  enh_filt <-
    d %>% 
    filter(baitID==enhID | preyID==enhID) %>%
    mutate(enhName=enhName)
  
  # filter out all the excluded delIDs
  for (delID in delIDs){
    enh_filt %<>% 
      filter(!baitID %in% delID &
               !preyID %in% delID)
  }

  
  # start the for enhancer loop
  #for (e in 1:nrow(enhData_new)){
    
    
  # get the IDs of all the enhancers
  enhIDs_all = data.frame(IDs=c(enhData_new$min_enhID,enhData_new$max_enhID)) %>% distinct() %>% arrange(IDs)
  
  # keep interactions between this enhancer and any other enhancer
  enh_filt %<>% 
    filter(baitID %in% enhIDs_all$IDs & preyID==enhID |
             preyID %in% enhIDs_all$IDs & baitID==enhID)
  # add to output df
  totals = rbind.data.frame(totals,enh_filt)

  
} # end of for enh loop




# get total E-E contacts in each group
# start by getting it for each sample
# and then split up the sample into group
sum_EE_contacts <-
  totals %>%
  group_by(sample,enhName) %>%
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
res.aov <- aov(contacts ~ tissue*genotype*promName*enhName, data = totals)
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
    res.aov <- aov(contacts ~ tissue*enh_numName, data = totals_temp)
    
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
      filter(term=="tissue:enh_numName" &
               adj.p.value<0.1) %>%
      separate(comparison, into=c("comp1","comp2"),sep="-") %>%
      separate(comp1, into=c("tissue1","enh1"),sep=":") %>%
      separate(comp2, into=c("tissue2","enh2"),sep=":") %>%
      filter(enh1==enh2)
  
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



####################################################
#### PLOT SEPARATE ENHANCERs ENHANCER CONTACTS #####
####################################################

# now calculate mean and sd of the sum of contacts for each group
# give the genotypes a shorter name
output <- 
  totals %>%
  group_by(tissue,genotype,enhName) %>%
  summarise(mean=mean(contacts),
            sd=sd(contacts)) %>%
  ungroup() %>%
  mutate(genotype=str_replace(genotype,"\\.CTCF\\.KO", "KO"))


# get the order of the cell types proper
output$tissue = factor(output$tissue, levels=c("Un", "Flk1", "CD41"))


# plot each enhancer separately
plt_list <- list()
for (i in 1:nrow(enhData_new)){ # for enh_numName loop
  
  # get the enhancer names
  enh=enhData_new$enh_numName[i]
  enhName=enhData_new$enhName[i]
  
  temp <-
  output %>%
    filter(enh_numName==enh)
  

  # plot these by cell type, and promoter
  plt_list[[enhName]] <-
    ggplot(temp, aes(x=tissue,y=mean))+
    geom_errorbar(aes(ymin=mean+sd, ymax=mean+sd),position=position_dodge(.9),width=.2)+
    #geom_linerange(aes(ymin = mean, ymax = mean+sd))+
    geom_bar(stat = "identity",position=position_dodge(), aes(fill=tissue))+
    scale_fill_brewer(palette="Accent")+
    #facet_grid(enhName~genotype)+
    ggtitle(enhName)+
    guides(fill=FALSE)+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid = element_blank(),
          strip.background =element_rect(fill="white"))
  
  
} # end of for enh_numName loop

  

plt_all <- plot_grid(plotlist=plt_list,ncol=5)


ggsave(plot=plt_all,
       filename=paste0(outFolder, "enhPromContacts_all.pdf"),
       width=12,
       height = 12)


#####################################
#### PLOT A SUBSET OF ENHANCERS #####
#####################################

# after looking at all the plots, find some to look at more closely
# enhancers to keep are:
# main list::: good_enh = c(-304,-58.4,-41.9,-43.1,3.2,24.2,47.5,63.6,64.6,106,110.2,204.7)

good_enh = c(-58.4,-41.9,3.2,24.2,110.2,204.7)


# filter the already calcualted mean contacts and sd in different cell type / genotype groups
enh_filt <-
  filter(output, enh_numName %in% good_enh)

# filter the enhancer data table as well
enhData_filt <-
  filter(enhData_new, enh_numName %in% good_enh)


# get the order of the cell types proper
enh_filt$tissue = factor(enh_filt$tissue, levels=c("Un", "Flk1", "CD41"))

# plot each enhancer separately
plt_list <- list()
for (i in 1:nrow(enhData_filt)){ # for enh_numName loop
  
  # get the enhancer names
  enhNum=enhData_filt$enh_numName[i]
  enhName=enhData_filt$enhName[i]
  
  temp <-
    enh_filt %>%
    filter(enh_numName==enhNum)
  
  for (i in 1:nrow(promData)){ # for prom loop
    
    # get the promoter names
    prom=strsplit(promData$promName[i],"\\.")[[1]][2]
    promID=promData$promID[i]
    
    enh_prom_filt <-
      temp %>%
      filter(promName==prom)
    
    # get a name for this plot combo
    combo = paste0(enhName," (",prom,")")
    
    # plot these by cell type, and promoter
    plt_list[[combo]] <-
      ggplot(enh_prom_filt, aes(x=tissue,y=mean))+
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
  
  
} # end of for enh_numName loop


plt_all <- plot_grid(plotlist=plt_list,ncol=4)


ggsave(plot=plt_all,
       filename=paste0(outFolder, "enhPromContacts_subset_smaller.pdf"),
       width=10,
       height = 3.5)

################################################
#### PLOT TOTAL ENHANCER PROMOTER CONTACTS #####
################################################

total_EPs <-
  totals %>%
  group_by(tissue,genotype) %>%
  summarise(mean=mean(contacts),
            sd=sd(contacts))

# get the order of the cell types proper
total_EPs$tissue = factor(total_EPs$tissue, levels=c("Un", "Flk1", "CD41"))


plot_EPs <-
  ggplot(total_EPs, aes(x=tissue,y=mean))+
  geom_errorbar(aes(ymin=mean+sd, ymax=mean+sd),position=position_dodge(.9),width=.2)+
  geom_linerange(aes(ymin = mean, ymax = mean+sd))+
  geom_bar(stat = "identity",position=position_dodge(), aes(fill=tissue))+
  scale_fill_brewer(palette="Accent")+
  facet_grid(~genotype)+
  ggtitle("Total Enhancer-Promoter Contacts")+
  #guides(fill=T)+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        strip.background =element_rect(fill="white"))


ggsave(plot=plot_EPs,
       filename=paste0(outFolder, "enhPromContacts_total_bothProms_sd.pdf"),
       width=6,
       height = 3)