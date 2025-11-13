# Compute score at the MLE

Compute score at the MLE

## Usage

``` r
score(mod, numerical = FALSE, get_score_covariance = FALSE)
```

## Arguments

- mod:

  an object of class `bbdml`

- numerical:

  Boolean. Defaults to `FALSE`. Indicator of whether to use the numeric
  Hessian and score (not recommended).

- get_score_covariance:

  Boolean. Defaults to `FALSE`. Should we return a robust estimate of
  variance of score: \\\hat{B}(\hat{\theta}) = \sum_i G(\hat{\theta};
  W_i) G(\hat{\theta}; W_i)^T\\. This parameter is not intended for
  users.

## Value

Score at the MLE. For \\G(\theta, w)\\ score function, returns \\\sum_i
G(\hat{\theta}, W_i)\\ if get_score_covariance = FALSE.

## Examples

``` r
data(soil_phylum_small_otu1)
mod <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)
score(mod)
#> [1]  6.339481e-08 -1.270095e-12  1.163048e-08 -5.889733e-12
```
