# mbkmeans: Mini-batch k-means clustering for single-cell RNA-seq

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![BioC release](http://www.bioconductor.org/shields/build/release/bioc/mbkmeans.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/mbkmeans)
[![BioC devel](http://www.bioconductor.org/shields/build/release/bioc/mbkmeans.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/mbkmeans)
[![R-CMD-check](https://github.com/drisso/mbkmeans/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/drisso/mbkmeans/actions)

This package implements the mini-batch k-means algorithm for large datasets, 
including support for on-disk data representation.

The method is described in details in the paper:

[S. Hicks, R. Liu, Y. Ni, E. Purdom, D. Risso (2021).
mbkmeans: Fast clustering for single cell data using mini-batch k-means. PLOS Computational Biology.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008625)

## Installation

In virtually all cases, installing from Bioconductor is recommended.

```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("mbkmeans")
```

In the rare event you need the development version from GitHub, use the following.

```{r}
library(devtools)
BiocManager::install("drisso/mbkmeans")
```
