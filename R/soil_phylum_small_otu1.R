#' Small soil phylum data for examples, sample data as data frame combined with counts for OTU 1 and sequencing depth.
#'
#' A small subset of \code{\link{soil_phylo_sample}} used for examples. A data frame made from the `phyloseq` object with only sample data and counts for OTU 1.
#'
#' @format A phyloseq-class experiment-level object with sample data and OTU 1 counts.
#' \describe{
#' \item{sam_data}{sample data with the following covariates:
#' \itemize{
#' \item \code{Plants}, values \code{0} and \code{1}. Index for different plants
#' \item \code{Day}, values \code{0} (initial sampling point), \code{1} (12 days after treatment additions), and \code{2} (82 days after treatment additions). Index for different days of measurement
#' \item \code{Amdmt}, values \code{0} (no additions), \code{1} (biochar additions), and \code{2} (fresh biomass additions). Index for different soil additives.
#' \item \code{DayAmdmt}, values \code{00}, \code{01}, \code{02}, \code{10}, \code{11}, \code{12}, \code{20}, \code{21}, and \code{22}. A single index for the combination of \code{Day} and \code{Amdmt} with \code{Day} as the first digit and \code{Amdmt} as the second digit.
#' \item \code{ID}, values \code{A}, \code{B}, \code{C}, \code{D},  and \code{F}. Index for different soil plots.
#' \item \code{W}, counts for OTU1 in each sample. This OTU corresponds with the phylum \emph{Proteobacteria}.
#' \item \code{M}, the sequencing depth for each sample.
#' }}
#' }
#' @references Whitman, T., Pepe-Ranney, C., Enders, A., Koechli, C., Campbell, A.,  Buckley, D. H., Lehmann, J. (2016). \emph{Dynamics of microbial community composi-tion and soil organic carbon mineralization in soil following addition of pyrogenic andfresh organic matter}. The ISME journal, 10(12):2918. <doi: 10.1038/ismej.2016.68>.
"soil_phylum_small_otu1"
