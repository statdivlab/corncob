# Identify differentially-abundant and differentially-variable taxa

Identify differentially-abundant and differentially-variable taxa

## Usage

``` r
differentialTest(
  formula,
  phi.formula,
  formula_null,
  phi.formula_null,
  data,
  link = "logit",
  phi.link = "logit",
  test,
  boot = FALSE,
  B = 1000,
  sample_data = NULL,
  taxa_are_rows = TRUE,
  filter_discriminant = TRUE,
  fdr_cutoff = 0.05,
  fdr = "fdr",
  full_output = FALSE,
  inits = NULL,
  inits_null = NULL,
  try_only = NULL,
  verbose = FALSE,
  robust = FALSE,
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

- formula_null:

  Formula for mean under null, without response

- phi.formula_null:

  Formula for overdispersion under null, without response

- data:

  a data frame containing the OTU table, `phyloseq`, or
  `SummarizedExperiment` object containing the variables in the models

- link:

  link function for abundance covariates, defaults to `"logit"`

- phi.link:

  link function for dispersion covariates, defaults to `"logit"`

- test:

  Character. Hypothesis testing procedure to use. One of `"Wald"`,
  `"LRT"` (likelihood ratio test), or `"Rao"`.

- boot:

  Boolean. Defaults to `FALSE`. Indicator of whether or not to use
  parametric bootstrap algorithm. (See
  [`pbWald`](https://statdivlab.github.io/corncob/reference/pbWald.md)
  and
  [`pbLRT`](https://statdivlab.github.io/corncob/reference/pbLRT.md)).

- B:

  Optional integer. Number of bootstrap iterations. Ignored if `boot` is
  `FALSE`. Otherwise, defaults to `1000`.

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

  Integer. Defaults to `0.05`. Desired false discovery rate.

- fdr:

  Character. Defaults to `"fdr"`. False discovery rate control method,
  see [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for more
  options.

- full_output:

  Boolean. Opetional. Defaults to `FALSE`. Indicator of whether to
  include full `bbdml` model output for all taxa.

- inits:

  Optional initializations for model fit using `formula` and
  `phi.formula` as rows of a matrix. Defaults to `NULL`.

- inits_null:

  Optional initializations for model fit using `formula_null` and
  `phi.formula_null` as rows of a matrix. Defaults to `NULL`.

- try_only:

  Optional numeric. Will try only the `try_only` taxa, specified either
  via numeric input or character taxa names. Useful for speed when
  troubleshooting. Defaults to `NULL`, testing all taxa.

- verbose:

  Boolean. Defaults to `FALSE`; print status updates for long-running
  analyses

- robust:

  Should robust standard errors be used? If not, model-based standard
  errors are used. Logical, defaults to `FALSE`.

- ...:

  Optional additional arguments for
  [`bbdml`](https://statdivlab.github.io/corncob/reference/bbdml.md)

## Value

An object of class `differentialTest`. List with elements `p` containing
the p-values, `p_fdr` containing the p-values after false discovery rate
control, `significant_taxa` containing the taxa names of the
statistically significant taxa, `significant_models` containing a list
of the model fits for the significant taxa, `all_models` containing a
list of the model fits for all taxa, `restrictions_DA` containing a list
of covariates that were tested for differential abundance,
`restrictions_DV` containing a list of covariates that were tested for
differential variability, `discriminant_taxa_DA` containing the taxa for
which at least one covariate associated with the abundance was perfectly
discriminant, `discriminant_taxa_DV` containing the taxa for which at
least one covariate associated with the dispersion was perfectly
discriminant, `data` containing the data used to fit the models. If
`full_output = TRUE`, it will also include `full_output`, a list of all
model output from `bbdml`.

## Details

See package vignette for details and example usage. Make sure the number
of columns in all of the initializations are correct! `inits` probably
shouldn't match `inits_null`. To use a contrast matrix, see
[`contrastsTest`](https://statdivlab.github.io/corncob/reference/contrastsTest.md).

## Examples

``` r
# data frame example
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

# phyloseq example (only run if you have phyloseq installed)
if (FALSE) { # \dontrun{
data_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylum_small_sample),
phyloseq::otu_table(soil_phylum_small_otu, taxa_are_rows = TRUE))
da_analysis <- differentialTest(formula = ~ DayAmdmt,
                               phi.formula = ~ DayAmdmt,
                               formula_null = ~ 1,
                               phi.formula_null = ~ DayAmdmt,
                               test = "Wald", boot = FALSE,
                               data = data_phylo,
                               fdr_cutoff = 0.05,
                               try_only = 1:5)
} # }
```
