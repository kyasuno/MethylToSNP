# R package: MethylToSNP

This is a modified version of [MethylToSNP](https://github.com/elnitskilab/MethylToSNP), which 
includes a major change in the mapping of indices of identified outliers of each cluster back to 
those in the original vector of cluster assignment. Additionally, for diagnostic purposes, 
plots of the distribution of betas are available when the data to be tested 
(probes with span >= 0.5 and the number of missing data `>= min.obs`) has <= 10 rows. 
Some parts of the code have also been vectorized. 
This version of package has been tested only in R 4.3.2. 

## Installation

To install this version:

```r
install.packages("devtools")
library(devtools)
install_github("kyasuno/MethylToSNP")
```
