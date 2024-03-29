#' differentialTest print function
#'
#' @param x Object of class \code{bbdml}
#' @param ... No optional arguments are accepted at this time.
#'
#'
#' @return \code{NULL}. Displays printed \code{differentialTest} summary.
#'
#' @examples
#' # phyloseq example
#' data(soil_phylum_small_sample)
#' data(soil_phylum_small_otu)
#' da_analysis <- differentialTest(formula = ~ DayAmdmt,
#'                                 phi.formula = ~ DayAmdmt,
#'                                 formula_null = ~ 1,
#'                                 phi.formula_null = ~ DayAmdmt,
#'                                 test = "Wald", boot = FALSE,
#'                                 data = soil_phylum_small_otu,
#'                                 sample_data = soil_phylum_small_sample,
#'                                 fdr_cutoff = 0.05,
#'                                 try_only = 1:5)
#' print(da_analysis)
#'
#' @export
print.differentialTest <- function(x, ...) {
  cat("Object of class differentialTest \n\n")
  cat("$p: p-values \n")
  cat("$p_fdr: FDR-adjusted p-values \n")
  cat("$significant_taxa: taxa names of the statistically significant taxa \n")
  cat("$significant_models: model summaries of the statistically significant taxa \n")
  cat("$all_models: all model summaries \n")
  cat("$restrictions_DA: covariates tested for differential abundance \n")
  cat("$restrictions_DV: covariates tested for differential variability \n")
  cat("$discriminant_taxa_DA: taxa for which at least one covariate associated with the abundance was perfectly discriminant \n")
  cat("$discriminant_taxa_DV: taxa for which at least one covariate associated with the dispersion was perfectly discriminant \n\n")
  cat("plot( ) to see a plot of tested coefficients from significant taxa \n")
}
