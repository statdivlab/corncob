#' Function to subset and convert phyloseq data
#'
#' @param data \code{phyloseq} object
#' @param select Name of OTU or taxa to select, must match taxa name in \code{data}
#'
#' @return A data.frame object
#' @examples
#' \dontrun{
#' data(soil_phylo)
#' convert_phylo(soil_phylo, "OTU1")
#' }
#'
#' @import phyloseq
#' @export
convert_phylo <- function(data, select) {
  subsamp <- phyloseq::prune_taxa(select, data)
  if (taxa_are_rows(input_data)) {
    W_tmp <- phyloseq::otu_table(subsamp)@.Data
  } else {
    W_tmp <- t(phyloseq::otu_table(subsamp)@.Data)
  }

  out <- data.frame(W = ,
                    M = phyloseq::sample_sums(data),
                    as.matrix(phyloseq::sample_data(subsamp)))
  return(out)
}
