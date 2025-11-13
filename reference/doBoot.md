# Function to run a bootstrap iteration

Internal function. Not intended for users.

## Usage

``` r
doBoot(mod, mod_null, test, robust = FALSE)
```

## Arguments

- mod:

  an object of class `bbdml`

- mod_null:

  an object of class `bbdml`

- test:

  Character. Hypothesis testing procedure to use. One of `"Wald"` or
  `"LRT"` (likelihood ratio test).

- robust:

  Should robust standard errors be used? If not, model-based standard
  arras are used. Logical, defaults to `FALSE`.

## Value

test statistic from one bootstrap iteration
