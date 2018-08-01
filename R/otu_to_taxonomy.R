#' Transform OTUs to their taxonomic label
#'
#' @param OTU String vector. Names of OTU labels in \code{data}
#' @param data \code{phyloseq} object with a taxonomy table
#'
#' @importFrom magrittr %>%
#'
#' @return String vector. Names of taxonomic labels matching labels of \code{OTU}.
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#'
#' @export
otu_to_taxonomy <- function(OTU, data) {
  return(apply(tax_table(data)[OTU,], 1, function(x) {paste(na.omit(x), collapse = '_')}))
}
