# ASPEN

A statistical method to analyze allele-specific expression (ASE) data. ASPEN estimates beta-binomial parameters, mean and dispersion, from the allelic counts.
ASPEN models allelic dispersion as a function of total gene expression and uses the Bayesian hierarchical model to shrink the estimated dispersion towards the expected one,
thereby stabilizing estimates and reducing the impact of spurious variability. 
ASPEN can be used to 
- identify genes with allelic imbalance (those with mean allelic ratio deviating from the theoretical value of 0.5);
- identify genes with allelic dispersion different from the expected one for the genes with similar expression levels;
- identify genes with changes in their allelic distribution over time;
- identify genes with changes in their allelic variation over time.

## Installation
```
if(!require(devtools)) install.packages("devtools")
library(devtools)

devtools::install_github("ewonglab/ASPEN")
```
