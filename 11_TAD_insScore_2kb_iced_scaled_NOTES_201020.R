





















## plot each TAD separately

for (tad in tadData$tad_name){
  
  # keep just this TAD
  insScores2 =
    insScores %>%
    filter(tad_name==!!tad)
  
  
  # get a variable for the ylimits
  lims=c((mean(insScores2$mean_ins)-.07),(mean(insScores2$mean_ins)+.07))
  
  
  ## change order of genotype variables
  insScores2$genotype=factor(insScores2$genotype, levels=c("WT", "P1.CTCF.KO", "P2.CTCF.KO"))
  insScores2$tissue=factor(insScores2$tissue, levels=c("Un", "Flk1", "CD41"))
  
  ## plot these insScores to see if there is a difference
  ggplot(insScores2, aes(x=group, y=mean_ins,fill=genotype))+
    stat_summary(geom="errorbar", fun.data=mean_se, width=.5)+
    stat_summary(geom="bar", fun="mean")+
    facet_grid(tad_name~tissue,scales="free")+
    scale_fill_brewer(palette = "Set1")+
    coord_cartesian(ylim = lims)+
    theme_bw()
  
  
  ## save the plot
  ggsave(filename=paste0(plotFolder, tad, "_mean_insScores_plot.pdf"),
         width=4,
         height=3)
  
  
} # end of for TAD plotting loop




## just WT plot for different cell type stages in fig2


for (tad in tadData$tad_name){
  
  # keep just this TAD and just WT samples
  insScores2 =
    insScores %>%
    filter(tad_name==!!tad & genotype=="WT")
  
  
  # get a variable for the ylimits
  lims=c((mean(insScores2$mean_ins)-.07),(mean(insScores2$mean_ins)+.07))
  
  
  ## change order of genotype variables
  insScores2$genotype=factor(insScores2$genotype, levels=c("WT", "P1.CTCF.KO", "P2.CTCF.KO"))
  insScores2$tissue=factor(insScores2$tissue, levels=c("Un", "Flk1", "CD41"))
  
  
  ## plot these insScores to see if there is a difference
  ggplot(insScores2, aes(x=tissue, y=mean_ins,fill=tissue))+
    stat_summary(geom="errorbar", fun.data=mean_se, width=.5)+
    stat_summary(geom="bar", fun="mean")+
    #facet_grid(tad_name,scales="free")+
    scale_fill_manual(values = c("#5f719dff", "#c5a0d2ff", "#599b78ff"))+
    coord_cartesian(ylim = lims)+
    theme_bw()
  
  
  ## save the plot
  ggsave(filename=paste0(plotFolder, tad, "_justWT_mean_insScores_plot.pdf"),
         width=4,
         height=3)
  
  
} # end of for TAD plotting loop





####### plot just main TAD alone ############



# keep just this TAD and just WT samples
insScores2 =
  insScores %>%
  filter(tad_name=="main_TAD" & genotype=="WT")


# get a variable for the ylimits
lims=c((mean(insScores2$mean_ins)-.07),(mean(insScores2$mean_ins)+.07))


## change order of genotype variables
insScores2$genotype=factor(insScores2$genotype, levels=c("WT", "P1.CTCF.KO", "P2.CTCF.KO"))
insScores2$tissue=factor(insScores2$tissue, levels=c("Un", "Flk1", "CD41"))


## plot these insScores to see if there is a difference
ggplot(insScores2, aes(x=tissue, y=mean_ins,fill=tissue))+
  stat_summary(geom="errorbar", fun.data=mean_se, width=.5)+
  stat_summary(geom="bar", fun="mean")+
  #facet_grid(tad_name,scales="free")+
  scale_fill_manual(values = c("#5f719dff", "#c5a0d2ff", "#599b78ff"))+
  coord_cartesian(ylim = lims)+
  theme_bw()


## save the plot
ggsave(filename=paste0(plotFolder, tad, "_justWT_mean_insScores_plot.pdf"),
       width=2.25,
       height=1.25)



####### plot just sub TADs alone ############



# keep just this TAD and just WT samples
insScores2 =
  insScores %>%
  filter((tad_name=="P1-P2_TAD" | tad_name=="P2-3'UTR_TAD")  & genotype=="WT")


# get a variable for the ylimits
lims=c((mean(insScores2$mean_ins)-.07),(mean(insScores2$mean_ins)+.07))


## change order of genotype variables
insScores2$genotype=factor(insScores2$genotype, levels=c("WT", "P1.CTCF.KO", "P2.CTCF.KO"))
insScores2$tissue=factor(insScores2$tissue, levels=c("Un", "Flk1", "CD41"))


## plot these insScores to see if there is a difference
ggplot(insScores2, aes(x=tissue, y=mean_ins,fill=tissue))+
  stat_summary(geom="errorbar", fun.data=mean_se, width=.5)+
  stat_summary(geom="bar", fun="mean")+
  facet_grid(~tad_name,scales="free")+
  scale_fill_manual(values = c("#5f719dff", "#c5a0d2ff", "#599b78ff"))+
  coord_cartesian(ylim = lims)+
  theme_bw()


## save the plot
ggsave(filename=paste0(plotFolder, tad, "_justWT_mean_insScores_plot.pdf"),
       width=2.8,
       height=1.25)





#### plot all TADs just for CD41 cells


for (tad in tadData$tad_name){
  
  # keep just this TAD and just WT samples
  insScores2 =
    insScores %>%
    filter(tad_name==!!tad & tissue=="CD41")
  
  
  # get a variable for the ylimits
  lims=c((mean(insScores2$mean_ins)-.04),(mean(insScores2$mean_ins)+.02))
  
  
  ## change order of genotype variables
  insScores2$genotype=factor(insScores2$genotype, levels=c("WT", "P1.CTCF.KO", "P2.CTCF.KO"))
  insScores2$tissue=factor(insScores2$tissue, levels=c("Un", "Flk1", "CD41"))
  
  
  ## plot these insScores to see if there is a difference
  ggplot(insScores2, aes(x=genotype, y=mean_ins,fill=genotype))+
    stat_summary(geom="errorbar", fun.data=mean_se, width=.5)+
    stat_summary(geom="bar", fun="mean")+
    #facet_grid(tad_name,scales="free")+
    scale_fill_brewer(palette = "Set1")+
    coord_cartesian(ylim = lims)+
    theme_bw()
  
  
  ## save the plot
  ggsave(filename=paste0(plotFolder, tad, "_justCD41_mean_insScores_plot.pdf"),
         width=3,
         height=1.5)
  
  
} # end of for TAD plotting loop




# calculate the mean and sd insulations scores by group
group_insScores =
  insScores %>%
  group_by(group) %>%
  summarise(mean_group_ins=mean(mean_ins),
            sd_group_ins=sd(mean_ins))



###########################################
#### ANOVA ON TOTAL PROMOTER CONTACTS #####
###########################################
## okay results so far, working code



# do stats (Three way anova)
#res.aov <- aov(mean_ins ~ tissue*genotype*tad_name, data = insScores)
#summary(res.aov)

## three way is not significant interaction effect tissue:genotype:tad_name

## now trying ONE WAY ANOVA just on WT samples and main TAD
### SOFIA SAYS THIS IS OKAY, AS THIS WAS THE ANALYSIS THAT WE WANTED IN THE FIRST PLACE
### SO ITS OKAY TO EXCLUDE THE OTHER GENOTYPES DATA AT THIS POINT
filt.dat = 
  insScores %>%
  filter(genotype=="WT" & tad_name=="main_TAD")


res.aov <- aov(mean_ins ~ tissue, data = filt.dat)
summary(res.aov)
tidy.res = tidy(res.aov)

# get the significant effects
sigEffects = tidy.res[tidy.res$p.value<0.05,]$term
sigEffects <- sigEffects[!is.na(sigEffects)]

# do the post-hocTest and keep only if that main / interaction effect is significant
phoc=
  tidy(TukeyHSD(res.aov)) %>%
  filter(term %in% sigEffects)



# export the ANOVA table
write.table(tidy.res,
            file=paste0(outFolder,"One-Way-ANOVA_mainTAD_WT.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            na = "",
            sep="\t")

# export the post Hoc results
write.table(phoc,
            file=paste0(outFolder,"One-Way-ANOVA_mainTAD_WT_postHoc.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            na = "",
            sep="\t")





### now doing a Two Way ANOVA for both sub TADs in WT cells only

### SO ITS OKAY TO EXCLUDE THE OTHER GENOTYPES DATA AT THIS POINT
filt.dat = 
  insScores %>%
  filter(genotype=="WT" & tad_name!="main_TAD")


res.aov <- aov(mean_ins ~ tissue*tad_name, data = filt.dat)
summary(res.aov)
tidy.res = tidy(res.aov)

# get the significant effects
sigEffects = tidy.res[tidy.res$p.value<0.05,]$term
sigEffects <- sigEffects[!is.na(sigEffects)]

# do the post-hocTest and keep only if that main / interaction effect is significant
phoc=
  tidy(TukeyHSD(res.aov)) %>%
  filter(term %in% sigEffects)





### fclean up the posthoc table and filter on interesting ones
phoc %<>%
  filter(term=="tissue:tad_name") %>%
  mutate(comparison=gsub("P1-P2","P1.P2",comparison)) %>%
  mutate(comparison=gsub("P2-3'UTR","P2.3'UTR",comparison)) %>%
  separate(comparison, sep="-", into=c("comp1","comp2")) %>%
  separate(comp1, sep=":", into=c("tissue1","TAD1")) %>%
  separate(comp2, sep=":", into=c("tissue2","TAD2")) %>%
  filter(TAD1==TAD2)





# export the ANOVA table
write.table(tidy.res,
            file=paste0(outFolder,"Two-Way-ANOVA_subTADs_WT.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            na = "",
            sep="\t")

# export the post Hoc results
write.table(phoc,
            file=paste0(outFolder,"Two-Way-ANOVA_subTADs_WT_postHoc.txt"),
            col.names=T,
            quote=F,
            row.names = F,
            na = "",
            sep="\t")







## now doing ANOVA on just CD41 cells all three TADs and genotypes

### SO ITS OKAY TO EXCLUDE THE OTHER GENOTYPES DATA AT THIS POINT
filt.dat = 
  insScores %>%
  filter(tissue=="CD41")


res.aov <- aov(mean_ins ~ genotype*tad_name, data = filt.dat)
summary(res.aov)
tidy.res = tidy(res.aov)


## not significant interaction effect p=0.719





## working up to here
#####################################################################################################



# TRIED THREE WAY ANOVA WITH CELL TYPE , GENOTYPE , AND TAD AS VARIABLES

## AND TRIED TWO WAY ANOVA WITH GROUP and TAD AS VARIABLES

### NON WERE SIGNIFICANT WITHIN CELL TYPE / TAD

# ie genotype does not have significant effect
