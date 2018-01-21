#' Overdispersion Gradient for a single sample
#'
#'
#' @keywords internal
#'
#' @export
gr_betast_i <- function(k, k_star, W, M, coth_st) {
  dg1 <- -digamma(M - W + (coth_st)/(k))
  dg2 <- -digamma((W + (k - 1)*(W + coth_st))/(k))
  dg3 <- digamma(((k - 1)*coth_st)/(k))
  dg4 <- digamma(M + coth_st)
  dg5 <- digamma(coth_st)
  dg6 <- digamma((coth_st)/(k))
  # dg7 <- digamma(M + coth_st) # This is dg4 again
  return( (4/(k*k_star*k_star)) * (k_star + 1) * (dg1 + (k - 1)*(dg2 + dg3 + dg4) - k*dg5 + dg6 + dg4))
}

#' Overdispersion Parameter Gradient Vector
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
#' @return Gradient of likelihood with respect to overdispersion parameters
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#'
#' @export
gr_betast <- function(theta, W, M, X, X_star, np, npstar, logpar = TRUE) {
  # extract matrix of betas (np x 1), first np entries
  b      <- utils::head(theta, np)
  # extract matrix of beta stars (npstar x 1), last npstar entries
  b_star <- utils::tail(theta, npstar)
  k      <- c(exp(X %*% b) + 1)
  k_star <- c(exp(2 * (X_star %*% b_star)) - 1)
  a2     <- 2 / (k * k_star)
  if (any(k_star == 0) || any(a2 < 0)) {
    # no overdispersion, density below has no b_star
    #val <- sum(dbinom(W, M, (k - 1)/k, log = TRUE))
    return(rep(0, npstar))
  }
  coth_st <- coth(X_star %*% b_star) - 1
  if (any(coth_st == 0)) {
    # no overdispersion, density below has no b_star
    #val <- sum(dbinom(W, M, (k - 1)/k, log = TRUE))
    out <- numDeriv::grad(func = dbetabin, x = theta,
                          W = W, M = M,
                          X = X,
                          X_star = X_star,
                          np = np,
                          npstar = npstar,
                          logpar = logpar)
    return(utils::tail(out, npstar))
  }

  if (!logpar) {
    stop("Use log.")
  }
  dg1 <- -digamma(M - W + (coth_st)/(k))
  dg2 <- -digamma((W + (k - 1)*(W + coth_st))/(k))
  dg3 <- digamma(((k - 1)*coth_st)/(k))
  dg4 <- digamma(M + coth_st)
  dg5 <- digamma(coth_st)
  dg6 <- digamma((coth_st)/(k))
  # dg7 <- digamma(M + coth_st) # This is dg4 again
  outNoX <- (4/(k*k_star*k_star)) * (k_star + 1) * (dg1 + (k - 1)*(dg2 + dg3 + dg4) - k*dg5 + dg6 + dg4)
  # outNoX <- mapply(gr_betast_i,
  #                  k = k, k_star = k_star, W = W, M = M,
  #                  coth_st = coth_st,
  #                  SIMPLIFY = TRUE)

  return(-c(crossprod(outNoX, X_star)))
}


