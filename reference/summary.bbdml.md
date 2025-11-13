# Summary function

Summary function

## Usage

``` r
# S3 method for class 'bbdml'
summary(object, ...)
```

## Arguments

- object:

  Object of class `bbdml`

- ...:

  No optional arguments are accepted at this time.

## Value

Object of class `summary.bbdml`. Displays printed model summary.

## Examples

``` r
data(soil_phylum_small_otu1)
mod <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)
summary(mod)
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
```
