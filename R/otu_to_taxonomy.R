#' Transform OTUs to their taxonomic label
#'
#' @param OTU String vector. Names of OTU labels in \code{data}
#' @param data \code{phyloseq} object with a taxonomy table
#' @param level (Optional). Character vector. Desired taxonomic levels for output.
#'
#' @importFrom magrittr %>%
#'
#' @return String vector. Names of taxonomic labels matching labels of \code{OTU}.
#'
#' @examples
#' \dontrun{
#' data(soil_phylo)
#' soil <- soil_phylo %>%
#' phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
#' phyloseq::tax_glom("Phylum")
#' da_analysis <- differentialTest(formula = ~ DayAmdmt,
#'                                 phi.formula = ~ DayAmdmt,
#'                                 formula_null = ~ 1,
#'                                 phi.formula_null = ~ DayAmdmt,
#'                                 test = "Wald", boot = FALSE,
#'                                 data = soil,
#'                                 fdr_cutoff = 0.05)
#' otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = soil, level = "Phylum")
#' }
#'
#' @export
otu_to_taxonomy <- function(OTU, data, level = NULL) {
  if (is.null(level)) {
    return(apply(tax_table(data)[OTU,], 1, function(x) {paste(stats::na.omit(x), collapse = '_')}))
  } else {
    return(apply(tax_table(data)[OTU, level], 1, function(x) {paste(stats::na.omit(x), collapse = '_')}))
  }
}
