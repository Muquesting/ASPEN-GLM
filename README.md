# ASPEN

A statistical method to analyze allele-specific expression (ASE) data. ASPEN estimates the mean and dispersion parameters of the beta-binomial distribution from allelic counts.
 Then, it models allelic dispersion as a function of total gene expression and shrinks the estimated dispersion toward the expected values using the Bayesian hierarchical framework. 
            
ASPEN can be used to 
- ASPEN can be used to (i) identify genes with allelic imbalance (those with mean allelic ratio deviating from the theoretical value of 0.5);
- identify genes with allelic dispersion deviating from the expected for the genes with similar expression levels;
- identify genes with changes in their allelic distribution across the groups;
- identify genes with changes in their allelic variation across the group.

## Installation
```
if(!require(devtools)) install.packages("devtools")
library(devtools)

devtools::install_github("ewonglab/ASPEN")
```
