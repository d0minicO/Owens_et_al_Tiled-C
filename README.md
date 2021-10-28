# Analysis of Tiled-C data
Instructions and custom scripts for anlaysis of Tiled-C data in our biorXiv paper here:

https://www.biorxiv.org/content/10.1101/2021.05.14.444178v1

any questions should be directed to dominic.owens@utoronto.ca


## Analysis part 1
The initial part of the analysis is based on Oudelaar et al 2020. Please see:

Oudelaar, A.M., Beagrie, R.A., Gosden, M. et al. Dynamics of the 4D genome during in vivo lineage specification and differentiation. Nat Commun 11, 2722 (2020). https://doi.org/10.1038/s41467-020-16598-7

This uses the Hi-CPro pipeline (Servant et al. Genome Biol 2015; https://github.com/nservant/HiC-Pro) 

Please see this link for additional details: https://github.com/oudelaar/TiledC

## Analysis part 2
*R code by Dominic Owens*:

**DISCLAIMER: The following R scripts are separated into sections with a description of what each one performs. The author developed this collection of scripts to meet his research team's specific data analysis needs and makes NO CLAIMS to their constituting a complete R package that will accept different inputs and reproducibly give a desired output across operating systems. These scripts are being shared in the spirit of open science to facilitate reproduction of the analysis done in our paper. These scripts could also serve as a useful starting point for user's own analysis. Inputs should be checked carefully as they may differ depending on your situation (eg number of samples, naming of samples, directory structures etc.) The outputs of each step MUST be checked to make sure they are being produced in a desired / expected format.**

#### Check and plot reporter numbers in each individual sample and merged samples

1_TileR_replicateAnalysis_ReporterNumbers_200322
1-2_TileR_2kb_raw_contactNumPerBaitPerSample_200404
1-3_TileR_2kb_raw_contactNumPerBaitPerSample_merged_200404


##### Performing clustering of individual samples and perform statistical testing (DESeq2) between genotypes and between cell types

 - 2_TileR_replicateAnalysis_DESeq_diffIntFrags_nonIced_matrices_v2_200324
 - 2-2_TileR_replicateAnalysis_DESeq_diffIntFrags_nonIced_matrices_justWT_200324
*(to plot PCA and heatmap of just WT samples)*


#### Extract normalised matrices for iced individual replicates to plot in python

 - 3_TileR_replicateAnalysis_iced_normMatrices_200327

#### Virtual Capture-C plots from viewpoints as specified in “baitsRightIDs.txt”

 - 4_TileR_2kb_iced_200331
*(bedGraph tracks outputted as hubs, including the replicate information as a SD (light colour above the track in UCSC))*

*Once individual replicates seem okay (clustering was okay and low sd on virtual CapC plots) then MERGE sam files into one sam file for each biological genotype/tissue like this*
`samtools merge mergedOutput.sam in1.sam in2.sam in3.sam` 

*Then follow the above Ice normalisation procedue as describe in* **Part 1** *to get a single ICE normalised matrix for each genotype/tissue*


#### Normalising individual replicate iced matrices and exporting matrices (whole region and smaller window)


#### Normalising merged replicate matrices and subtracting to allow comparison between genotypes and cell types


#### Filtering matrices to only include statistically significantly different bins calculated by DESeq, and subsetting matrices to only smaller window

