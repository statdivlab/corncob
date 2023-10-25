#' Get quantiles of beta binom
#'
#' @param p Numeric. Probability for quantile
#' @param M Numeric vector of sequencing depth
#' @param mu Numeric vector of abundance parameter
#' @param phi Numeric vector of dispersion parameter
#'
#' @return quantile
#'
#' @examples
#' data(soil_phylum_small)
#' mod <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small)
#' qbetabinom(.5, M = mod$M[1], mu = mod$mu.resp[1], phi = mod$phi.resp[1])
#' @export
qbetabinom <- function(p, M, mu, phi)  {
  pdfvec <- dbetabinom_cts(x = 0:M, size = M, prob = mu, rho = phi)
  # get cumulative probability
  CDF <- cumsum(pdfvec)
  # Find first value above quantile cutoff
  above <- which(CDF > p)[1]
  # Subtract one
  quant <- above - 1

  if (quant > 0) {

    return(quant + (p - CDF[quant])/pdfvec[above])
  } else {
    return(p/pdfvec[1])
  }
}
