# Parametric bootstrap Wald test

Parametric bootstrap Wald test

## Usage

``` r
pbWald(mod, mod_null, B = 1000, robust = FALSE)
```

## Arguments

- mod:

  an object of class `bbdml`

- mod_null:

  an object of class `bbdml`, should be nested within `mod`

- B:

  Integer. Defaults to `1000`. Number of bootstrap iterations.

- robust:

  Should robust standard errors be used? If not, model-based standard
  arras are used. Logical, defaults to `FALSE`.

## Value

P-value from parametric bootstrap Wald test.

## Examples

``` r
data(soil_phylum_small_otu1)
mod1 <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)

mod2 <- bbdml(formula = cbind(W, M - W) ~ 1,
phi.formula = ~ 1,
data = soil_phylum_small_otu1)
pbWald(mod1, mod2, B = 50)
#> [1] 0.01960784
```
