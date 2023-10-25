#' Negative betabinomial density
#'
#' Created as a convenient helper function for optimization. Not intended for users.
#'
#' @param theta Numeric vector. Parameters associated with \code{X} and \code{X_star}
#' @param W Numeric vector of counts
#' @param M Numeric vector of sequencing depth
#' @param X Matrix of covariates associated with abundance (including intercept)
#' @param X_star Matrix of covariates associated with dispersion (including intercept)
#' @param np Number of covariates associated with abundance (including intercept)
#' @param npstar Number of covariates associated with dispersion (including intercept)
#' @param link ink function for abundance covariates
#' @param phi.link ink function for dispersion covariates
#' @param logpar Boolean. Defaults to \code{TRUE}. Indicator of whether to return log-likelihood.
#'
#' @return Negative beta-binomial (log-)likelihood
dbetabin_neg <- function(theta, W, M, X, X_star, np, npstar, link, phi.link, logpar = TRUE) {

  # extract matrix of betas (np x 1), first np entries
  b      <- utils::head(theta, np)
  # extract matrix of beta stars (npstar x 1), last npstar entries
  b_star <- utils::tail(theta, npstar)

  mu.withlink <- X %*% b
  phi.withlink <- X_star %*% b_star
  mu <- switch(link, "logit" = invlogit(mu.withlink))
  phi <- switch(phi.link, "fishZ" = invfishZ(phi.withlink), "logit" = invlogit(phi.withlink))

  val <- sum(dbetabinom_cts(W, M, prob = mu, rho = phi, log = TRUE))

  return(-val)
}


#' Betabinomial density
#'
#' @param theta Numeric vector. Parameters associated with \code{X} and \code{X_star}
#' @param W Numeric vector of counts
#' @param M Numeric vector of sequencing depth
#' @param X Matrix of covariates associated with abundance (including intercept)
#' @param X_star Matrix of covariates associated with dispersion (including intercept)
#' @param np Number of covariates associated with abundance (including intercept)
#' @param npstar Number of covariates associated with dispersion (including intercept)
#' @param link ink function for abundance covariates
#' @param phi.link ink function for dispersion covariates
#' @param logpar Boolean. Defaults to \code{TRUE}. Indicator of whether to return log-likelihood.
#'
#' @return Negative beta-binomial (log-)likelihood
dbetabin <- function(theta, W, M, X, X_star, np, npstar, link, phi.link, logpar = TRUE) {

  # extract matrix of betas (np x 1), first np entries
  b      <- utils::head(theta, np)
  # extract matrix of beta stars (npstar x 1), last npstar entries
  b_star <- utils::tail(theta, npstar)

  mu.withlink <- X %*% b
  phi.withlink <- X_star %*% b_star
  mu <- switch(link, "logit" = invlogit(mu.withlink))
  phi <- switch(phi.link, "fishZ" = invfishZ(phi.withlink), "logit" = invlogit(phi.withlink))

  val <- sum(dbetabinom_cts(W, M, prob = mu, rho = phi, log = TRUE))

  return(val)
}
