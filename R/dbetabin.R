#' Betabinomial density for a single sample
#'
#'
#' @keywords internal
#'
#' @export
dbetabin_i <- function(a1, a2, W, M, logpar = TRUE) {
  if (logpar) {
    return(lbeta(a1 + W, a2 + M - W) - lbeta(a1, a2) + lchoose(M, W))
  }
}

#' Betabinomial density
#'
#' @param theta parameters
#' @param W absolute abundance
#' @param M sample size
#' @param X mean covariates
#' @param X_star overdispersion covariates
#' @param np number of mean parameters
#' @param npstar number of overdisperion parameters
#' @param logpar Indicator of log-likelihood, defaults to TRUE
#'
#' @return Beta-binomial log-likelihood
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#'
#' @export
dbetabin <- function(theta, W, M, X, X_star, np, npstar, logpar = TRUE) {

  # extract matrix of betas (np x 1), first np entries
  b <- utils::head(theta, np)


  # extract matrix of beta stars (npstar x 1), last npstar entries
  b_star <- utils::tail(theta, npstar)

  k      <- c(exp(X %*% b) + 1)
  k_star <- c(exp(2 * (X_star %*% b_star)) - 1)
  a2     <- 2 / (k * k_star)
  if (sum(a2) == Inf || any(a2 < 0)) {
    # no overdispersion
    val <- sum(dbinom(W, M, (k - 1)/k, log = TRUE))
    return(-val)
  }
  a1     <- a2 * (k - 1)

  # val <- sum(mapply(dbetabin_i,
  #                   a1 = a1, a2 = a2, W = W, M = M,
  #                   MoreArgs = list(logpar = logpar)))
  val <- sum(lbeta(a1 + W, a2 + M - W) - lbeta(a1, a2) + lchoose(M, W))

  return(-val)
}
