# Get highest density interval of beta-binomial

Get highest density interval of beta-binomial

## Usage

``` r
HDIbetabinom(percent, M, mu, phi)
```

## Arguments

- percent:

  Numeric. Percent interval desired.

- M:

  Numeric vector of sequencing depth

- mu:

  Numeric vector of abundance parameter

- phi:

  Numeric vector of dispersion parameter

## Value

List where `lower` represents the lower bound and `upper` represents the
upper bound

## Examples

``` r
data(soil_phylum_small_otu1)
mod <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)
HDIbetabinom(.95, M = mod$M[1], mu = mod$mu.resp[1], phi = mod$phi.resp[1])
#> $lower
#> [1] 61719
#> 
#> $upper
#> [1] 74743
#> 
```
