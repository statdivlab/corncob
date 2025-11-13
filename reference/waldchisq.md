# Wald-type chi-squared test

Wald-type chi-squared test

## Usage

``` r
waldchisq(
  mod,
  mod_null = NULL,
  restrictions = NULL,
  restrictions.phi = NULL,
  contrasts_DA = NULL,
  contrasts_DV = NULL,
  robust = FALSE
)
```

## Arguments

- mod:

  an object of class `bbdml`

- mod_null:

  Optional. An object of class `bbdml`, should be nested within `mod`.
  If not included, need to include `restrictions` or `restrictions.phi`.

- restrictions:

  Optional. Defaults to `NULL`. Numeric vector indicating the parameters
  associated with the abundance to test, or character vector with name
  of variable to test. Note that `1` is the intercept associated with
  the abundance.

- restrictions.phi:

  Optional. Defaults to `NULL`. Numeric vector indicating the parameters
  associated with the dispersion to test, or character vector with name
  of variable to test. Note that `1` is the intercept associated with
  the dispersion.

- contrasts_DA:

  List. Optional. Constructs a contrast matrix. List elements should be
  characters specifying contrasts in the parameters within `formula`.
  Note that this is only available with `"Wald"` value for `test`.

- contrasts_DV:

  List. Optional. Constructs a contrast matrix. List elements should be
  characters specifying contrasts in the parameters within
  `phi.formula`. Note that this is only available with `"Wald"` value
  for `test`.

- robust:

  Should robust standard errors be used? If not, model-based standard
  arras are used. Logical, defaults to `FALSE`.

## Value

Matrix with wald test statistics and p-values. Only performs univariate
tests.

P-value from Wald test.

## Examples

``` r
data(soil_phylum_small_otu1)
mod1 <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)

mod2 <- bbdml(formula = cbind(W, M - W) ~ 1,
phi.formula = ~ 1,
data = soil_phylum_small_otu1)

# Example using mod_null
waldchisq(mod = mod1, mod_null = mod2)
#> [1] 6.60754e-06
#> attr(,"df")
#> [1] 2
waldchisq(mod = mod1, mod_null = mod2, robust = TRUE)
#> [1] 0.0001905122
#> attr(,"df")
#> [1] 2

# Example using restrictions and restrictions.phi
waldchisq(mod = mod1, restrictions = 2, restrictions.phi = 2)
#> [1] 6.60754e-06
#> attr(,"df")
#> [1] 2
waldchisq(mod = mod1, restrictions = "DayAmdmt", restrictions.phi = "DayAmdmt")
#> [1] 6.60754e-06
#> attr(,"df")
#> [1] 2
waldchisq(mod = mod1, restrictions = 2, restrictions.phi = "DayAmdmt")
#> [1] 6.60754e-06
#> attr(,"df")
#> [1] 2
waldchisq(mod = mod1, restrictions = 2, restrictions.phi = 2, robust = TRUE)
#> [1] 0.0001905122
#> attr(,"df")
#> [1] 2
```
