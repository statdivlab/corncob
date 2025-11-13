# Maximum Likelihood for the Beta-binomial Distribution

Maximum Likelihood for the Beta-binomial Distribution

## Usage

``` r
bbdml(
  formula,
  phi.formula,
  data,
  link = "logit",
  phi.link = "logit",
  method = "trust",
  control = list(maxit = 1000, reltol = 1e-14),
  numerical = FALSE,
  nstart = 1,
  inits = NULL,
  allow_noninteger = FALSE,
  robust = FALSE,
  ...
)
```

## Arguments

- formula:

  an object of class `formula`: a symbolic description of the model to
  be fitted to the abundance

- phi.formula:

  an object of class `formula` without the response: a symbolic
  description of the model to be fitted to the dispersion

- data:

  a data frame, `phyloseq`, or `SummarizedExperiment` object containing
  the variables in the models

- link:

  link function for abundance covariates, defaults to `"logit"`

- phi.link:

  link function for dispersion covariates, defaults to `"logit"`

- method:

  optimization method, defaults to `"trust"`, or see
  [`optimr`](https://rdrr.io/pkg/optimx/man/optimr.html) for other
  options

- control:

  optimization control parameters (see
  [`optimr`](https://rdrr.io/pkg/optimx/man/optimr.html))

- numerical:

  Boolean. Defaults to `FALSE`. Indicator of whether to use the numeric
  Hessian (not recommended).

- nstart:

  Integer. Defaults to `1`. Number of starts for optimization.

- inits:

  Optional initializations as rows of a matrix. Defaults to `NULL`.

- allow_noninteger:

  Boolean. Defaults to `FALSE`. Should noninteger W's and M's be
  allowed? This behavior was not permitted prior to v4.1, needs to be
  explicitly allowed.

- robust:

  Should robust standard errors be returned? If not, model-based
  standard arras are used. Logical, defaults to `FALSE`.

- ...:

  Optional additional arguments for
  [`optimr`](https://rdrr.io/pkg/optimx/man/optimr.html) or
  [`trust`](https://rdrr.io/pkg/trust/man/trust.html)

## Value

An object of class `bbdml`.

## Examples

``` r
# data frame example
data(soil_phylum_small_otu1)
bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)
#> 
#> Call:
#> bbdml(formula = cbind(W, M - W) ~ DayAmdmt, phi.formula = ~DayAmdmt, 
#>     data = soil_phylum_small_otu1)
#> 
#> 
#> Coefficients associated with abundance:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -0.44595    0.03604 -12.375 7.18e-13 ***
#> DayAmdmt21  -0.16791    0.04067  -4.129 0.000297 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> 
#> Coefficients associated with dispersion:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  -5.3077     0.3537 -15.008 6.44e-15 ***
#> DayAmdmt21   -1.3518     0.5029  -2.688    0.012 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> 
#> Log-likelihood: -286.53

# phyloseq example (only run this if you have phyloseq installed)
if (FALSE) { # \dontrun{
data(soil_phylum_small_sample)
data(soil_phylum_small_otu)
data_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylum_small_sample),
phyloseq::otu_table(soil_phylum_small_otu, taxa_are_rows = TRUE))
bbdml(formula = Proteobacteria ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = data_phylo)
} # }
```
