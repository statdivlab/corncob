# Parametric bootstrap likelihood ratio test

Parametric bootstrap likelihood ratio test

## Usage

``` r
pbLRT(mod, mod_null, B = 1000)
```

## Arguments

- mod:

  an object of class `bbdml`

- mod_null:

  an object of class `bbdml`, should be nested within `mod`

- B:

  Integer. Defaults to `1000`. Number of bootstrap iterations.

## Value

P-value from parametric bootstrap likelihood ratio test.

## Examples

``` r
data(soil_phylum_small_otu1)
mod1 <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)

mod2 <- bbdml(formula = cbind(W, M - W) ~ 1,
phi.formula = ~ 1,
data = soil_phylum_small_otu1)
pbLRT(mod1, mod2, B = 50)
#> [1] 0.01960784
```
