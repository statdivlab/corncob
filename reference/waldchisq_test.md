# Wald-type chi-squared test statistic (model-based or robust)

This is a helper function and not intended for users

## Usage

``` r
waldchisq_test(
  mod,
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

Test statistic for Wald test.
