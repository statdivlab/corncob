# differentialTest print function

differentialTest print function

## Usage

``` r
# S3 method for class 'differentialTest'
print(x, ...)
```

## Arguments

- x:

  Object of class `bbdml`

- ...:

  No optional arguments are accepted at this time.

## Value

`NULL`. Displays printed `differentialTest` summary.

## Examples

``` r
# phyloseq example
data(soil_phylum_small_sample)
data(soil_phylum_small_otu)
da_analysis <- differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ DayAmdmt,
                                test = "Wald", boot = FALSE,
                                data = soil_phylum_small_otu,
                                sample_data = soil_phylum_small_sample,
                                fdr_cutoff = 0.05,
                                try_only = 1:5)
print(da_analysis)
#> Object of class differentialTest 
#> 
#> $p: p-values 
#> $p_fdr: FDR-adjusted p-values 
#> $significant_taxa: taxa names of the statistically significant taxa 
#> $significant_models: model summaries of the statistically significant taxa 
#> $all_models: all model summaries 
#> $restrictions_DA: covariates tested for differential abundance 
#> $restrictions_DV: covariates tested for differential variability 
#> $discriminant_taxa_DA: taxa for which at least one covariate associated with the abundance was perfectly discriminant 
#> $discriminant_taxa_DV: taxa for which at least one covariate associated with the dispersion was perfectly discriminant 
#> 
#> plot( ) to see a plot of tested coefficients from significant taxa 
```
