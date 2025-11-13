# Get index of restricted terms for Wald test

Created as a convenient helper function. Not intended for users.

## Usage

``` r
getRestrictionTerms(
  mod,
  mod_null = NULL,
  restrictions = NULL,
  restrictions.phi = NULL
)
```

## Arguments

- mod:

  an object of class `bbdml`

- mod_null:

  Optional. An object of class `bbdml`. Defaults to `NULL`

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

## Value

A list with `mu` representing the index of the restricted covariates
associated with abundance and `phi` representing the index of the
restricted covarates associated with the dispersion
