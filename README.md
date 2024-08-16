# Collecting and Identification the Outbreak Cluster

<!-- badges: start -->

<!-- badges: end -->

The goal of Collecting and Identification the Outbreak Cluster `caIRA` is to find the cluster based tree and metadata of the data. The package is based on the paper:

Ragonnet-Cronin, M., Hodcroft, E., Hué, S. et al. Automated analysis of phylogenetic clusters. BMC Bioinformatics 14, 317 (2013). <https://doi.org/10.1186/1471-2105-14-317>

Hall M, Woolhouse M, Rambaut A (2015) Epidemic Reconstruction in a Phylogenetics Framework: Transmission Trees as Partitions of the Node Set. PLoS Comput Biol 11(12): e1004613. <https://doi.org/10.1371/journal.pcbi.1004613>

This package is the part of Dhihram Tenrisau, MSc Health Data Science summer project, 'Phylodynamic of Norovirus in UK 2003-2023'. The project is supervised by Stéphane Hué

``` r
# install.packages("devtools")
devtools::install_github("Dhihram/huebreaker")
```

his package needs the additional package `tidyverse`, `ape`, `treeio`, and `dplyr`

``` r
library(tidyverse)
library(ape)
library(dplyr)
library(treeio)
library(caIRA)
```

You need the 2 data in this file, the first data is tree files (newick or nexus) and metadata file. The metadata file consists of `label`, `location`, and `date` columns. The `label` column must be same with the label in tree file.

This package with `genclus` will utilize:
1. Finding and clustering the monophylectic groups in the tree
2. Add the parameter of the clusters: `bootstrap_treshold`, `data_range`, and `samearea`
3. Keep the maximum monophylectic groups in the cluster identify

the `bootstrap_treshold` is the minimum bootstrap value to be considered as a cluster. The `data_range` is the range of the days to be considered as a cluster. The `samearea` is the boolean value to consider the same area as a cluster.

``` r
res <- genclus(tree, metat, bootstrap_threshold = 80, date_range = 30, samearea = TRUE)
```
For the manual, you can see [here](https://dhihram.github.io/caIRA/)
