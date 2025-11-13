# Compute sandwich standard errors. Legacy function. Use sand_vcov instead.

Compute sandwich standard errors. Legacy function. Use sand_vcov
instead.

## Usage

``` r
sandSE(mod, numerical = FALSE)
```

## Arguments

- mod:

  an object of class `bbdml`

- numerical:

  Boolean. Defaults to `FALSE`. Indicator of whether to use the numeric
  Hessian and score (not recommended).

## Value

Sandwich variance-covariance matrix

## Examples

``` r
data(soil_phylum_small_otu1)
mod <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)
sandSE(mod)
#> Warning: You called sandSE(). Please use sand_vcov() instead. They are the same function but sand_vcov() has more transparent naming.
#>              (Intercept)   DayAmdmt21  (Intercept)   DayAmdmt21
#> (Intercept)  0.001318049 -0.001318049  0.006663799 -0.006663799
#> DayAmdmt21  -0.001318049  0.001678019 -0.006663799  0.010663491
#> (Intercept)  0.006663799 -0.006663799  0.156697831 -0.156697831
#> DayAmdmt21  -0.006663799  0.010663491 -0.156697831  0.313975615
```
