# Compute Hessian matrix at the MLE

Compute Hessian matrix at the MLE

## Usage

``` r
hessian(mod, numerical = FALSE)
```

## Arguments

- mod:

  an object of class `bbdml`

- numerical:

  Boolean. Defaults to `FALSE`. Indicator of whether to use the numeric
  Hessian (not recommended).

## Value

Hessian matrix at the MLE. In this setting, it's hard to compute
expectations to calculate the information matrix, so we return the
consistent estimate using sample moments: \\\hat{A}(\hat{\theta}) =
\sum_i \frac{\partial^2}{\partial \theta \partial \theta^T} l(\theta,
W_i)\\ evaluated at \\\theta = \hat{\theta}\\.

## Examples

``` r
data(soil_phylum_small_otu1)
mod <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)
hessian(mod)
#>             (Intercept)  DayAmdmt21 (Intercept) DayAmdmt21
#> (Intercept) 3587.083848 2816.656102   -4.254135  -2.521716
#> DayAmdmt21  2816.656102 2816.656102   -2.521716  -2.521716
#> (Intercept)   -4.254135   -2.521716   15.826817   7.827544
#> DayAmdmt21    -2.521716   -2.521716    7.827544   7.827544
```
