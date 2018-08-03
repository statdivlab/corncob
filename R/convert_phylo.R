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
  W_tmp <- matrix(phyloseq::otu_table(subsamp)@.Data, ncol = 1)


  out <- data.frame(W = W_tmp,
                    M = phyloseq::sample_sums(data),
                    data.frame(phyloseq::sample_data(subsamp)))
  return(out)
}
