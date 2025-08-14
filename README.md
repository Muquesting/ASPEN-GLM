# ASPEN

A statistical method to analyze allele-specific expression (ASE) data. ASPEN estimates the mean and dispersion parameters of the beta-binomial distribution from allelic counts.
 Then, it models allelic dispersion as a function of total gene expression and shrinks the estimated dispersion toward the expected values using the Bayesian hierarchical framework. 
            
ASPEN can be used to 
- identify genes with allelic imbalance;
- identify genes with allelic dispersion deviating from the expected for the genes with similar expression levels;
- identify genes with changes in their allelic distribution across the groups;
- identify genes with changes in their allelic variation across the group.

## Installation
```
if(!require(devtools)) install.packages("devtools")
library(devtools)

devtools::install_github("ewonglab/ASPEN")
```

## Quick start
```
library(ASPEN)
#loading reference allele and total count matrices
data("Bl6_Cast_a1")
data("Bl6_Cast_tot")

#estimating beta-binomial distribution parameters
bb_init_params <- estim_bbparams(Cast_B6_a1, Cast_B6_tot, min_cells = 5, cores = 6)

#shrinking dispersion
shrink_pars <- correct_theta(bb_init_params, delta_set = 50, N_set = 30, thetaFilter = 0.001)

#preparing a list of genes that are excluded from the reference allele mapping bias evaluation
load_file <- system.file("extdata", "mm10_genesXY.txt", package = "ASPEN")
genesXY <- read.table(load_file)
load_file <- system.file("extdata", "mm10_imprinted_genes.xlsx", package = "ASPEN")
genesIMPR <- read.xlsx(load_file, colNames = T)
genes2remove <- c(genesXY$V1, genesIMPR$imprinted.genes)

#evaluating reference allele bias
global_estims <- glob_disp(Cast_B6_a1, Cast_B6_tot, genes.excl = genes2remove, min_counts = 5)

bb_mean_res <- bb_mean(Cast_B6_a1, Cast_B6_tot,
                       shrink_pars, min_cells = 5,
                       min_counts = 5, glob_params = global_estims)
```

## Tutorials

For the step-by-step guide on how to use ASPEN, please see [ASPEN tutorial](https://ewonglab.github.io/ASPEN/) 

## Code from the main manuscript

A repository with the code to reproduce analyses from the ASPEN [manuscript](https://www.biorxiv.org/content/10.1101/2025.04.16.649227v1) can be found [here](https://github.com/ewonglab/ASPEN_manuscript/)

