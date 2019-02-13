#' Rename taxa
#'
#' Renames taxa to have short human-readable names
#'
#' @param x Object of class \code{phyloseq}
#' @param name Character, defaults to \code{"OTU"}. Optional. String to use in every taxa name.
#'
#' @details The original taxa names are saved as the \code{original_names} attribute. See the example for an example of how to access the original names.
#'
#' @return Object of class \code{phyloseq}, with taxa renamed (defaults to OTU1, OTU2, ...), with the original taxa names saved as an attribute.
#'
#' @examples
#' \dontrun{
#' data(soil_phylo)
#' x <- clean_taxa_names(soil_phylo)
#' # Use this line to see the original taxa names
#' attr(x, "original_names")
#' }
#' @export
clean_taxa_names <- function(x, name = "OTU") {
  if ("phyloseq" %in% class(x)) {
    attr(x, "original_names") <- phyloseq::taxa_names(x)
    phyloseq::taxa_names(x) <- paste0(name, seq(phyloseq::ntaxa(x)))
    return(x)
  } else {
    stop("clean_taxa_names is intended for phyloseq objects!")
  }
}
