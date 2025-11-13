# Parameter Gradient Vector

Used for internal optimization. Not intended for users.

## Usage

``` r
gr_full(theta, W, M, X, X_star, np, npstar, link, phi.link, logpar = TRUE)
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

Gradient of likelihood with respect to parameters
