#' Get highest density interval of beta-binomial
#'
#' @param percent Numeric. Percent interval desired.
#' @param M Numeric vector of sequencing depth
#' @param mu Numeric vector of abundance parameter
#' @param phi Numeric vector of dispersion parameter
#'
#' @return List where \code{lower} represents the lower bound and \code{upper} represents the upper bound
#'
#' @examples
#' data(soil_phylum_small)
#' mod <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small)
#' HDIbetabinom(.95, M = mod$M[1], mu = mod$mu.resp[1], phi = mod$phi.resp[1])
#'
#' @export
HDIbetabinom <- function(percent, M, mu, phi)  {
  pdfvec <- VGAM::dbetabinom(x = 0:M, size = M, prob = mu, rho = phi)

  myord <- order(pdfvec, decreasing = TRUE)
  # Get ordered values
  ordered_pdf <- pdfvec[myord]
  # get cumulative probability
  ordered_CDF <- cumsum(ordered_pdf)
  # Find first value above ordered cutoff
  above <- which(ordered_CDF > percent)[1]

  numbers_in_HDI <- c(0:M)[myord][1:above]

  return(list(lower = min(numbers_in_HDI), upper = max(numbers_in_HDI)))
}
