# Analysis of Tiled-C data

[![DOI](https://zenodo.org/badge/393514046.svg)](https://zenodo.org/badge/latestdoi/393514046)

Instructions and custom scripts for anlaysis of Tiled-C data in our Nature Communications paper here:

https://www.nature.com/articles/s41467-022-28376-8

Details are available in the methods section of the paper. Any questions that arise should be directed to dominic.owens@utoronto.ca


## Analysis part 1
The initial part of the analysis is based on Oudelaar et al 2020. Please see:

Oudelaar, A.M., Beagrie, R.A., Gosden, M. et al. Dynamics of the 4D genome during in vivo lineage specification and differentiation. Nat Commun 11, 2722 (2020). https://doi.org/10.1038/s41467-020-16598-7

This uses the Hi-CPro pipeline (Servant et al. Genome Biol 2015; https://github.com/nservant/HiC-Pro) 

Please see this link for additional details: https://github.com/oudelaar/TiledC

## Analysis part 2

Raw matrix data and quality control plots and statistics are included in this github submission. Raw sequencing data has been uploaded to GEO. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184490

*R code by Dominic Owens*:

**DISCLAIMER: The following R scripts are separated into sections with a description of what each one performs. The author developed this collection of scripts to meet his research team's specific data analysis needs and makes NO CLAIMS to their constituting a complete R package that will accept different inputs and reproducibly give a desired output across operating systems. These scripts are being shared in the spirit of open science and to facilitate reproduction of the analyses done in our paper. These scripts could also serve as a useful starting point for user's own analysis. Inputs should be checked carefully as they may differ depending on your situation (eg number of samples, naming of samples, directory structures etc.) The outputs of each step MUST be checked to make sure they are being produced in a desired / expected format.**

#### 1. Check and plot reporter numbers in each individual sample and merged samples

 - 1_TileR_replicateAnalysis_ReporterNumbers_200322
 - 1-2_TileR_2kb_raw_contactNumPerBaitPerSample_200404
 - 1-3_TileR_2kb_raw_contactNumPerBaitPerSample_merged_200404 *(to determine and plot reporter numbers PER BAIT in each individual sample and merged samples)*


#### 2. Performing clustering of individual samples and perform statistical testing (DESeq2) between genotypes and between cell types

 - 2_TileR_replicateAnalysis_DESeq_diffIntFrags_nonIced_matrices_v2_200324
 - 2-2_TileR_replicateAnalysis_DESeq_diffIntFrags_nonIced_matrices_justWT_200324
*(to plot PCA and heatmap of just WT samples)*


#### 3. Extract normalised matrices for iced individual replicates to plot in python

 - 3_TileR_replicateAnalysis_iced_normMatrices_200327

#### 4. Virtual Capture-C plots from viewpoints as specified in ???baitsRightIDs.txt???

 - 4_TileR_2kb_iced_200331
*(bedGraph tracks outputted as hubs, including the replicate information as a SD (light colour above the track in UCSC))*

*Once individual replicates seem okay (clustering was okay and low sd on virtual CapC plots) then MERGE sam files into one sam file for each biological genotype/tissue like this*

```samtools merge mergedOutput.sam in1.sam in2.sam in3.sam```

*Then follow the above ICE normalisation procedue as describe in* **Part 1** *to get a single ICE normalised matrix for each genotype/tissue*


#### 5. Calculate and plot reporter numbers in merged samples

 - 5_TileR_replicateAnalysis_ReporterNumbers_merged_libs_200322 *(Input is merged and ICE normalised matrices)*


#### 6. Normalising merged replicate matrices and subtracting to allow comparison between genotypes and cell types

 - 6_TileR_subtractionAndSubsetR_DESeq_filt_200322

*Input is merged iced matrices*
*Outputs are*
-- normalised matrices over whole tiled window
-- normalised matrices over smaller specified window
-- subtracted matrices over whole tiled window
-- subtracted matrices over smaller specified window
-- DESeq2 differential bins over whole tiled window
-- DESeq2 differential bins over smaller specified window

#### 7. Calculate total contacts from P1 and P2 promoters 

 - 7_totalPromContacts_2kb_iced_scaled_200422

*Input is iced and scaled matrices*


#### 8. Calculate contacts from P1 and P2 promoters to enhancers

 - 8_enhPromContacts_2kb_iced_scaled_200423


#### 9. Calculate contacts from P1 and P2 promoters to CTCF sites

 - 9_ctcfPromContacts_2kb_iced_scaled_200424


#### 10. Calculate TAD insulation scores

 - 10_TAD_insScore_2kb_iced_scaled_201020

#### 11. Calculate contacts between CTCF sites

 - 11_ctcfCTCFContacts_2kb_iced_scaled_201117




## Analysis part 3 -- Visualisation

Plotting matrices code was orignally written by Marieke Oudelaar. See the file here: https://github.com/oudelaar/TiledC/blob/master/TiledC_matrix_visualisation.py

Further development was done to automate the process for multiple matrices in a directory and to allow custom colour maps:

 - mainPlotter_relThresh_customCmap

*Make sure to set your directories accordingly and to set the n_bins to the right number. You can derive this from the coordinates.bed file generated by Tiled_sam2rawmatrix_MO.pl. You should also specify a threshold for plotting (95th is a good starting point).*
