# Rao-type chi-squared test (model-based or robust)

Rao-type chi-squared test (model-based or robust)

## Usage

``` r
raotest(mod, mod_null)
```

## Arguments

- mod:

  an object of class `bbdml`

- mod_null:

  an object of class `bbdml`, should be nested within `mod`

## Value

P-value from likelihood ratio test.

## Examples

``` r
data(soil_phylum_small_otu1)
mod1 <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)

mod2 <- bbdml(formula = cbind(W, M - W) ~ 1,
phi.formula = ~ 1,
data = soil_phylum_small_otu1)
raotest(mod1, mod2)
#> [1] 1
```
