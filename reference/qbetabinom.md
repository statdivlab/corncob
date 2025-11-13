# Get quantiles of beta binom

Get quantiles of beta binom

## Usage

``` r
qbetabinom(p, M, mu, phi)
```

## Arguments

- p:

  Numeric. Probability for quantile

- M:

  Numeric vector of sequencing depth

- mu:

  Numeric vector of abundance parameter

- phi:

  Numeric vector of dispersion parameter

## Value

quantile

## Examples

``` r
data(soil_phylum_small_otu1)
mod <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)
qbetabinom(.5, M = mod$M[1], mu = mod$mu.resp[1], phi = mod$phi.resp[1])
#> [1] 68185.71
```
