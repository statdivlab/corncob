# Generate initialization for optimization

Generate initialization for optimization

## Usage

``` r
genInits(W, M, X, X_star, np, npstar, link, phi.link, nstart = 1, use = TRUE)
```

## Arguments

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

- nstart:

  Integer. Defaults to `1`. Number of starts for optimization.

- use:

  Boolean. Defaults to `TRUE`. Indicator of whether to use deterministic
  intialization.

## Value

Matrix of initializations

## Examples

``` r
set.seed(1)
seq_depth <- rpois(20, lambda = 10000)
my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
my_covariate <- cbind(rep(c(0,1), each = 10))
colnames(my_covariate) <- c("X1")

genInits(W = my_counts, M = seq_depth,
       X = cbind(1, my_covariate), X_star = cbind(1, my_covariate),
       np = 2, npstar = 2,
       link = "logit",
       phi.link = "logit", nstart = 2, use = TRUE)
#>                        X1                    
#> [1,] -4.669910 0.03427153  0.0000000 0.000000
#> [2,] -4.381932 0.28231162 -0.1073783 0.309576
```
