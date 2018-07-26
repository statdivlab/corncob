#' Get quantiles of beta binom
#'
#' @param p probability for quantile
#' @param size size
#' @param mu probability parameter
#' @param phi overdispersion parameter
#'
#' @return quantile
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
qbetabinom <- function(p, size, mu, phi)  {
  pdfvec <- VGAM::dbetabinom(x = 0:size, size = size, prob = mu, rho = phi)
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
