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
#' @export
otu_to_taxonomy <- function(OTU, data, level = NULL) {
  if (requireNamespace("phyloseq", quietly = TRUE)) {
    if ("phyloseq" %in% class(data)) {
      if (is.null(level)) {
        return(apply(phyloseq::tax_table(data)[OTU,], 1, function(x) {paste(stats::na.omit(x), collapse = '_')}))
      } else {
        return(apply(phyloseq::tax_table(data)[OTU, level], 1, function(x) {paste(stats::na.omit(x), collapse = '_')}))
      }
    } else {
      stop("This function currently only works for phyloseq objects.")
    }
  } else {
    warn_phyloseq()
  }
}


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
#' @export
otu_to_taxonomy_se <- function(OTU, data, level = NULL) {
  if (!inherits(data, "SummarizedExperiment")) {
    stop("This function currently only works for SummarizedExperiment objects.")
  }

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    warn_sumexp()
  }

  if (is.null(level)) {
    apply(SummarizedExperiment::rowData(data)[OTU,], 1, function(x) {paste(stats::na.omit(x), collapse = "_")})
  } else {
    sapply(SummarizedExperiment::rowData(data)[OTU, level], function(x) {paste(stats::na.omit(x), collapse = "_")})
  }
}
