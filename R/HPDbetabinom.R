#' Get HPD interval of beta binom
#'
#' @param percent percent interval
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
HPDbetabinom <- function(percent, size, mu, phi)  {
  pdfvec <- VGAM::dbetabinom(x = 0:size, size = size, prob = mu, rho = phi)

  myord <- order(pdfvec, decreasing = TRUE)
  # Get ordered values
  ordered_pdf <- pdfvec[myord]
  # get cumulative probability
  ordered_CDF <- cumsum(ordered_pdf)
  # Find first value above ordered cutoff
  above <- which(ordered_CDF > percent)[1]

  numbers_in_HPD <- c(0:size)[myord][1:above]

  return(list(lower = min(numbers_in_HPD), upper = max(numbers_in_HPD)))
}
