#' Function to subset and convert phyloseq data
#'
#' @param data a \code{phyloseq} object
#' @param select Name of OTU or taxa to select, must match taxa name in \code{data}
#'
#' @return A \code{data.frame} object, with elements \code{W} as the observed counts, \code{M} as the sequencing depth, and the sample data with their original names.
#'
#' @export
convert_phylo <- function(data, select) {
  if (requireNamespace("phyloseq", quietly = TRUE)) {
    subsamp <- suppressWarnings(phyloseq::prune_taxa(select, data))
    W_tmp <- matrix(phyloseq::otu_table(subsamp)@.Data, ncol = 1)


    out <- data.frame(W = W_tmp,
                      M = phyloseq::sample_sums(data),
                      data.frame(phyloseq::sample_data(subsamp)))
    return(out)
  } else {
    warn_phyloseq()
  }
}
