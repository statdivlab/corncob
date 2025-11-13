# Objective function

Used for internal optimization. Not intended for users.

## Usage

``` r
objfun(theta, W, M, X, X_star, np, npstar, link, phi.link)
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

## Value

List of negative log-likelihood, gradient, and hessian
