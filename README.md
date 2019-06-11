# corncob <img src="docs/logo.png" align="right" width="165px"/>
Count Regression for Correlated Observations with the Beta-binomial

[![Build Status](https://travis-ci.org/bryandmartin/corncob.svg?branch=master)](https://travis-ci.org/bryandmartin/corncob)
[![codecov](https://codecov.io/gh/bryandmartin/CORNCOB/branch/master/graph/badge.svg?token=GnLFG7QNsh)](https://codecov.io/gh/bryandmartin/CORNCOB)

## Installation

To download the corncob package, use the code below.

``` r
# install.packages("devtools")
devtools::install_github("bryandmartin/corncob")
library(corncob)
```

## Use

The vignette demonstrates example usage of all main functions. You can see the vignette by using the following code:

``` r
# install.packages("devtools")
devtools::install_github("bryandmartin/corncob", build_vignettes = TRUE)
library(corncob)
# Use this to view the vignette in the corncob HTML help
help(package = "corncob", help_type = "html")
# Use this to view the vignette as an isolated HTML file
utils::browseVignettes(package = "corncob")
```
## Status

The preprint describing the methodology implemented in corncob is available [here](https://arxiv.org/abs/1902.02776).

## Bug Reports / Change Requests

If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/bryandmartin/corncob/issues).
