#' Function to throw error if the `phyloseq` package is called but it is not installed
#'
#' @export
warn_phyloseq <- function() {
  stop("You are trying to use a `phyloseq` data object or `phyloseq` helper function without having the `phyloseq` package installed. Please either install the package or use a standard data frame.")
}
