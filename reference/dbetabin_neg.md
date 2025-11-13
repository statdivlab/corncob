# Negative betabinomial density

Created as a convenient helper function for optimization. Not intended
for users.

## Usage

``` r
dbetabin_neg(theta, W, M, X, X_star, np, npstar, link, phi.link, logpar = TRUE)
```

## Arguments

- theta:

  Numeric vector. Parameters associated with `X` and `X_star`

- W:

  Numeric vector of counts

- M:

  Numeric vector of sequencing depth

- X:

  Matrix of covariates associated with abundance (including intercept)

- X_star:

  Matrix of covariates associated with dispersion (including intercept)

- np:

  Number of covariates associated with abundance (including intercept)

- npstar:

  Number of covariates associated with dispersion (including intercept)

- link:

  ink function for abundance covariates

- phi.link:

  ink function for dispersion covariates

- logpar:

  Boolean. Defaults to `TRUE`. Indicator of whether to return
  log-likelihood.

## Value

Negative beta-binomial (log-)likelihood
