# Check for nested models

Check for nested models

## Usage

``` r
checkNested(mod, mod_null)
```

## Arguments

- mod:

  an object of class `bbdml`

- mod_null:

  an object of class `bbdml`

## Value

`TRUE` if `mod_null` is nested within `mod`, otherwise it throws an
error.

## Examples

``` r
data(soil_phylum_small_otu1)
mod1 <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)

mod2 <- bbdml(formula = cbind(W, M - W) ~ 1,
phi.formula = ~ 1,
data = soil_phylum_small_otu1)

checkNested(mod1, mod2)
#> [1] TRUE
```
