# Densities of beta binomial distributions, permitting non integer x and size

In some cases we may not have integer W and M's. In these cases, we can
still use corncob to estimate parameters, but we need to think of them
as no longer coming from the specific beta binomial parametric model,
and instead from an estimating equations framework.

## Usage

``` r
dbetabinom_cts(x, size, prob, rho = 0, log = FALSE)
```

## Arguments

- x:

  the value at which defined the density

- size:

  number of trials

- prob:

  the probability of success

- rho:

  the correlation parameter

- log:

  if TRUE, log-densities p are given

## Author

Thomas W Yee

Xiangjie Xue

Amy D Willis
