# Identify differentially-abundant and differentially-variable taxa using contrasts

Identify differentially-abundant and differentially-variable taxa using
contrasts

## Usage

``` r
contrastsTest(
  formula,
  phi.formula,
  contrasts_DA = NULL,
  contrasts_DV = NULL,
  data,
  link = "logit",
  phi.link = "logit",
  sample_data = NULL,
  taxa_are_rows = TRUE,
  filter_discriminant = TRUE,
  fdr_cutoff = 0.05,
  fdr = "fdr",
  inits = NULL,
  try_only = NULL,
  ...
)
```

## Arguments

- formula:

  an object of class `formula` without the response: a symbolic
  description of the model to be fitted to the abundance

- phi.formula:

  an object of class `formula` without the response: a symbolic
  description of the model to be fitted to the dispersion

- contrasts_DA:

  List. Optional. Constructs a contrast matrix. List elements should be
  characters specifying contrasts in the parameters within `formula`.
  Note that this is only available with `"Wald"` value for `test`. Must
  include at least one of `contrasts_DA` or `contrasts_DV`.

- contrasts_DV:

  List. Optional. Constructs a contrast matrix. List elements should be
  characters specifying contrasts in the parameters within
  `phi.formula`. Note that this is only available with `"Wald"` value
  for `test`. Must include at least one of `contrasts_DA` or
  `contrasts_DV`.

- data:

  a data frame containing the OTU table, or `phyloseq` object containing
  the variables in the models

- link:

  link function for abundance covariates, defaults to `"logit"`

- phi.link:

  link function for dispersion covariates, defaults to `"logit"`

- sample_data:

  Data frame or matrix. Defaults to `NULL`. If `data` is a data frame or
  matrix, this must be included as covariates/sample data.

- taxa_are_rows:

  Boolean. Optional. If `data` is a data frame or matrix, this indicates
  whether taxa are rows. Defaults to `TRUE`.

- filter_discriminant:

  Boolean. Defaults to `TRUE`. If `FALSE`, discriminant taxa will not be
  filtered out.

- fdr_cutoff:

  Integer. Defaults to `0.05`. Desired type 1 error rate

- fdr:

  Character. Defaults to `"fdr"`. False discovery rate control method,
  see [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for more
  options.

- inits:

  Optional initializations for model fit using `formula` and
  `phi.formula` as rows of a matrix. Defaults to `NULL`.

- try_only:

  Optional numeric. Will try only the `try_only` taxa, specified either
  via numeric input or character taxa names. Useful for speed when
  troubleshooting. Defaults to `NULL`, testing all taxa.

- ...:

  Optional additional arguments for
  [`bbdml`](https://statdivlab.github.io/corncob/reference/bbdml.md)

## Value

An object of class `contrastsTest`. List with elements `p` containing
the p-values for each contrast, `p_fdr` containing the p-values after
false discovery rate control, `significant_taxa` containing the taxa
names of the statistically significant taxa, `contrasts_DA` containing
the contrast matrix for parameters associated with the abundance,
`contrasts_DV` containing the contrast matrix for parameters associated
with the dispersion, `discriminant_taxa_DA` containing the taxa for
which at least one covariate associated with the abundance was perfectly
discriminant, `discriminant_taxa_DV` containing the taxa for which at
least one covariate associated with the dispersion was perfectly
discriminant, and `data` containing the data used to fit the models.

## Details

This function uses contrast matrices to test for differential abundance
and differential variability using a Wald-type chi-squared test. To use
a formula implementation, see
[`differentialTest`](https://statdivlab.github.io/corncob/reference/differentialTest.md).

## Examples

``` r
# data frame example
# note that this function will only run if the `limma` package is installed
limma_install <- try(find.package("limma"), silent = TRUE)
if (!(inherits(limma_install, "try-error"))) {
  data(soil_phylum_contrasts_sample)
  data(soil_phylum_contrasts_otu)
  da_analysis <- contrastsTest(formula = ~ DayAmdmt,
                              phi.formula = ~ DayAmdmt,
                              contrasts_DA = list("DayAmdmt21 - DayAmdmt11",
                                                   "DayAmdmt22 - DayAmdmt21"),
                              data = soil_phylum_contrasts_otu,
                               sample_data = soil_phylum_contrasts_sample,
                               fdr_cutoff = 0.05,
                              try_only = 1:5)
}
#> Renaming (Intercept) to Intercept
#> Renaming (Intercept) to Intercept
#> Renaming (Intercept) to Intercept
#> Renaming (Intercept) to Intercept
#> Renaming (Intercept) to Intercept
#> Renaming (Intercept) to Intercept
#> Renaming (Intercept) to Intercept
#> Renaming (Intercept) to Intercept
#> Renaming (Intercept) to Intercept
#> Renaming (Intercept) to Intercept

# phyloseq example (only run if you have phyloseq installed)
if (FALSE) { # \dontrun{
contrasts_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylum_contrasts_sample),
phyloseq::otu_table(soil_phylum_contrasts_otu, taxa_are_rows = TRUE))
da_analysis <- contrastsTest(formula = ~ DayAmdmt,
                             phi.formula = ~ DayAmdmt,
                             contrasts_DA = list("DayAmdmt21 - DayAmdmt11",
                                                 "DayAmdmt22 - DayAmdmt21"),
                             data = contrasts_phylo,
                             fdr_cutoff = 0.05,
                             try_only = 1:5)
} # }
```
