# Parametric bootstrap Rao test

Parametric bootstrap Rao test

## Usage

``` r
pbRao(mod, mod_null, B = 1000)
```

## Arguments

- mod:

  an object of class `bbdml`

- mod_null:

  an object of class `bbdml`, should be nested within `mod`

- B:

  Integer. Defaults to `1000`. Number of bootstrap iterations.

## Value

P-value from parametric bootstrap Rao test.

## Examples

``` r
data(soil_phylum_small_otu1)
mod1 <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)

mod2 <- bbdml(formula = cbind(W, M - W) ~ 1,
phi.formula = ~ 1,
data = soil_phylum_small_otu1)
pbRao(mod1, mod2, B = 10)
#> [1] 1
```
