#' Rename taxa
#'
#' Renames taxa to have short human-readable names
#'
#' @param x Object of class \code{phyloseq}
#' @param name Character, defaults to \code{"OTU"}. Optional. String to use in every taxa name.
#'
#' @return Object of class \code{phyloseq}, with taxa renamed (defaults to OTU1, OTU2, ...)
#' @export
clean_taxa_names <- function(x, name = "OTU") {
  if ("phyloseq" %in% class(x)) {
    taxa_names(x) <- paste0(name, seq(phyloseq::ntaxa(x)))
    return(x)
  } else {
    stop("clean_taxa_names is intended for phyloseq objects!")
  }
}
