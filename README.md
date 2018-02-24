
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
--------

Capturing uncertainty in genetic population assignment is one of the key advantages of using Bayesian inference through programs such as [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html). However, this uncertainty is often ignored when using assignment inferences in downstream analysis, such as in landscape genetics.

The `STUCTURE.popgen` package overcomes these limitations by providing functions for using supplied q-matrix assignment probabilities to account for uncertainty while calculating population genetics metrics among inferred clusters. Functions are also provided for working more seamlessly through well-established `R` packages (`adegenet`, `PopGenReport`, and `hierfstat`) to calculate individual-based and population-based population genetics metrics.

Installation
------------

You can install `STUCTURE.popgen` from `GitHub` with:

``` r
# install.packages("devtools")
devtools::install_github("nicholasjclark/STRUCTURE.popgen")
```

Brief Introduction to Current Functions
---------------------------------------

The `STRUCTURE.popgen` function is the workhorse of the package. This function randomly assigns individuals to pre-identified population clusters using STRUCTURE assignment probabilities, which are provided in a supplied q-matrix (Pritchard et al. 2000).

A range of wrapper functions are also provided to calculate various commonly-used population genetics statistics. For example, the `fis.STRUCTURE.popgen` function assigns individuals to clusters using STRUCTURE assignment probabilities provided in the qmatrix, calculates individual *Fis* values and compares cluster-level *Fis* values using an analysis of variance (ANOVA). This can be repeated many times (using the `nreps` argument) to calculate summary statistics (i.e. confidence intervals) by allowing for uncertainty in individual assignment probabilities. Other functions are currently provided to calculate *Fst*, allelic richness and inbreeding probabilities.

This package is currently under development, so watch this space to keep track of improved functionality and flexibility.

References
----------

Pritchard, J.K., Stephens, M. & Donnelly, P. (2000) Inference of population structure using multilocus genotype data. *Genetics*, 155, 945-959.
