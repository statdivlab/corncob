# Compute sandwich estimate of variance-covariance matrix

Compute sandwich estimate of variance-covariance matrix

## Usage

``` r
sand_vcov(mod, numerical = FALSE)
```

## Arguments

- mod:

  an object of class `bbdml`

- numerical:

  Boolean. Defaults to `FALSE`. Indicator of whether to use the numeric
  Hessian and score (not recommended).

## Value

Sandwich variance-covariance matrix. \\\hat{A}^{-1} \hat{B}
\hat{A}^{-1}\\.

## Examples

``` r
data(soil_phylum_small_otu1)
mod <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)
sand_vcov(mod)
#>              (Intercept)   DayAmdmt21  (Intercept)   DayAmdmt21
#> (Intercept)  0.001318049 -0.001318049  0.006663799 -0.006663799
#> DayAmdmt21  -0.001318049  0.001678019 -0.006663799  0.010663491
#> (Intercept)  0.006663799 -0.006663799  0.156697831 -0.156697831
#> DayAmdmt21  -0.006663799  0.010663491 -0.156697831  0.313975615
```
