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

Documentation is now complete and includes examples. We are currently working on expanding vignettes and tutorials. Please leave an issue if you have any requests for a tutorial.

Manuscript is in preparation.

## Bug Reports / Change Requests

If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/bryandmartin/corncob/issues).
