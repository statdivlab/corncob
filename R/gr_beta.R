#' Beta Gradient for a single sample
#'
#'
#' @keywords internal
#'
#' @export
gr_beta_i <- function(k, k_star, W, M, coth_st) {
  dg1 <- -digamma(M - W + coth_st/(k))
  dg2 <- digamma((W + (k - 1)*(W + coth_st))/(k))
  dg3 <- digamma((coth_st)/(k))
  dg4 <- -digamma(((k - 1)*coth_st)/(k))
  return( (1/(k^2)) * (k - 1) * coth_st * (dg1 + dg2 + dg3 + dg4))
}

#' Mean Parameter Gradient Vector
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
#' @return Gradient of likelihood with respect to mean parameters
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#'
#' @export
gr_beta <- function(theta, W, M, X, X_star, np, npstar, logpar = TRUE) {
  # extract matrix of betas (np x 1), first np entries
  b      <- utils::head(theta, np)
  # extract matrix of beta stars (npstar x 1), last npstar entries
  b_star <- utils::tail(theta, npstar)
  k      <- c(exp(X %*% b) + 1)
  k_star <- c(exp(2 * (X_star %*% b_star)) - 1)

  coth_st <- coth(X_star %*% b_star) - 1
  if (!logpar) {
    stop("Use log.")
  }
  outNoX <- mapply(gr_beta_i,
                   k = k, k_star = k_star, W = W, M = M,
                   coth_st = coth_st,
                   SIMPLIFY = TRUE)
  # Will be n-vector, want out np vector
  # Should be np gradients, one for each beta
  # X is n by np, need colsums after multiplying by row. Cross product
  return(-c(crossprod(outNoX, X)))
}


