#' Function to subset and convert SummarizedExperiment data
#'
#' @param data a \code{SummarizedExperiment} object
#' @param select Name of OTU or taxa to select, must match taxa name in \code{data}
#'
#' @return A \code{data.frame} object, with elements \code{W} as the observed counts, \code{M} as the sequencing depth, and the sample data with their original names.
#'
#' @export
convert_sumexp <- function(data, select) {
  if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    subsamp <- data[select, ]
    W_tmp <- matrix(t(SummarizedExperiment::assay(subsamp)), ncol = 1)


    out <- data.frame(W = W_tmp,
                      M = colSums(SummarizedExperiment::assay(data)),
                      SummarizedExperiment::colData(subsamp))
    return(out)
  } else {
    warn_sumexp()
  }
}
