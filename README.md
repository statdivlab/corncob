# corncob <img src="docs/logo.png" align="right" width="165px"/>
Count Regression for Correlated Observations with the Beta-binomial

`corncob` is an `R` package for modeling relative abundance and testing hypotheses about the effect of covariates on relative abundance. The `corncob` methodology was specifically developed for modelling microbial abundances based on high throughput sequencing data, such as 16S or whole-genome sequencing.

[![Build Status](https://travis-ci.com/bryandmartin/corncob.svg?branch=master)](https://travis-ci.com/bryandmartin/corncob)
[![codecov](https://codecov.io/gh/bryandmartin/CORNCOB/branch/master/graph/badge.svg?token=GnLFG7QNsh)](https://codecov.io/gh/bryandmartin/CORNCOB)
[![Docker Repository on Quay](https://quay.io/repository/fhcrc-microbiome/corncob/status "Docker Repository on Quay")](https://quay.io/repository/fhcrc-microbiome/corncob)

## Installation

To download the corncob package, use the code below.

``` r
# install.packages("devtools")
devtools::install_github("bryandmartin/corncob")
library(corncob)
```

## Docker

Instead of installing corncob to your local system, you can use corncob via the pre-compiled Docker image: `quay.io/fhcrc-microbiome/corncob`.


## Use

The vignette demonstrates example usage of all main functions. Please [file an issue](https://github.com/bryandmartin/corncob/issues) if you have a request for a tutorial that is not currently included. You can see the vignette by using the following code (note that this requires a TeX installation to view properly):


``` r
# install.packages("devtools")
devtools::install_github("bryandmartin/corncob", build_vignette = TRUE, build_opts = c())
library(corncob)
# Use this to view the vignette in the corncob HTML help
help(package = "corncob", help_type = "html")
# Use this to view the vignette as an isolated HTML file
utils::browseVignettes(package = "corncob")
```

Note that R does not allow variable names to start with numbers. Sometimes, when going directly from QIIME2 to phyloseq objects, taxa names will be a large string starting with numbers. To clean these taxa names for use with corncob, use  `clean_taxa_names(my_phyloseq_object)`, see `?clean_taxa_names` for details.

## Citation

If you use `corncob` for your analysis, please cite our manuscript:

Bryan D. Martin, Daniela Witten, and Amy D. Willis. (2020). [*Modeling microbial abundances and dysbiosis with beta-binomial regression*.](https://projecteuclid.org/euclid.aoas/1587002666) Annals of Applied Statistics, Volume 14, Number 1, pages 94-115.

An open-access preprint is available on arXiv [here](https://arxiv.org/abs/1902.02776).

## Bug Reports / Change Requests

If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/bryandmartin/corncob/issues).
