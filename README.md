# corncob <img src="docs/logo.png" align="right" width="165px"/>
Count Regression for Correlated Observations with the Beta-binomial

[![Build Status](https://travis-ci.org/bryandmartin/CORNCOB.svg?branch=master)](https://travis-ci.org/bryandmartin/CORNCOB)
[![codecov](https://codecov.io/gh/bryandmartin/CORNCOB/branch/master/graph/badge.svg?token=GnLFG7QNsh)](https://codecov.io/gh/bryandmartin/CORNCOB)

## Installation

To download the corncob package, use the code below.

``` r
# install.packages("devtools")
devtools::install_github("bryandmartin/corncob")
library(corncob)
```

## Use

There is currently a vignette created as a tutorial for the STAMPS 2018 workshop. 
This vignette demonstrates example usage of all main functions. You can see the vignette by using the following code:

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

Documentation is complete, but not finalized, and most examples are still in progress. However, you can see examples of most functions intended for users in the package vignette.

Manuscript is in preparation.

## Bug Reports / Change Requests

If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/bryandmartin/corncob/issues).
