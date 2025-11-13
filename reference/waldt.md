# Wald-type t test (model-based or robust)

Wald-type t test (model-based or robust)

## Usage

``` r
waldt(mod)
```

## Arguments

- mod:

  an object of class `bbdml`

## Value

Matrix with wald test statistics and p-values. Only performs univariate
tests.

## Examples

``` r
data(soil_phylum_small_otu1)
mod <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)
waldt(mod)
#>                   Estimate Std. Error    t value     Pr(>|t|)
#> mu.(Intercept)  -0.4459452 0.03603627 -12.374899 7.184737e-13
#> mu.DayAmdmt21   -0.1679129 0.04066628  -4.129045 2.970223e-04
#> phi.(Intercept) -5.3077030 0.35365558 -15.008113 6.439048e-15
#> phi.DayAmdmt21  -1.3517509 0.50285497  -2.688153 1.195984e-02
```
