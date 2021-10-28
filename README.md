# Analysis of Tiled-C data
Instructions and custom scripts for anlaysis of Tiled-C data in our biorXiv paper here:

https://www.biorxiv.org/content/10.1101/2021.05.14.444178v1


## Analysis part 1
The initial part of the analysis is based on Oudelaar et al 2020. Please see:

Oudelaar, A.M., Beagrie, R.A., Gosden, M. et al. Dynamics of the 4D genome during in vivo lineage specification and differentiation. Nat Commun 11, 2722 (2020). https://doi.org/10.1038/s41467-020-16598-7

This uses the Hi-CPro pipeline (Servant et al. Genome Biol 2015; https://github.com/nservant/HiC-Pro) 

Please see this link for additional details: https://github.com/oudelaar/TiledC

## Analysis part 2
*R code by Dominic Owens*: 
- Check and plot reporter numbers in each individual sample and merged samples
- Performing clustering of individual samples (DESeq2)
- Statistical testing between genotypes and between cell types
- Normalising individual replicate iced matrices and exporting matrices (whole region and smaller window)
- Normalising merged replicate matrices and subtracting to allow comparison between genotypes and cell types
- Filtering matrices to only include statistically significantly different bins calculated by DESeq, and subsetting matrices to only smaller window

