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

  val <- sum(VGAM::dbetabinom(W, M, prob = mu, rho = phi, log = TRUE))

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
#'
#' @examples
#' set.seed(1)
#' seq_depth <- rpois(20, lambda = 10000)
#' my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
#' my_covariate <- cbind(rep(c(0,1), each = 10))
#' colnames(my_covariate) <- c("X1")
#'
#' example_data <- data.frame("W" = my_counts, "M" = seq_depth, my_covariate)
#'
#' dbetabin(theta = rep(-4, 4), W = my_counts, M = seq_depth,
#'        X = cbind(1, my_covariate), X_star = cbind(1, my_covariate),
#'        np = 2, npstar = 2,
#'        link = "logit",
#'        phi.link = "logit")
dbetabin <- function(theta, W, M, X, X_star, np, npstar, link, phi.link, logpar = TRUE) {

  # extract matrix of betas (np x 1), first np entries
  b      <- utils::head(theta, np)
  # extract matrix of beta stars (npstar x 1), last npstar entries
  b_star <- utils::tail(theta, npstar)

  mu.withlink <- X %*% b
  phi.withlink <- X_star %*% b_star
  mu <- switch(link, "logit" = invlogit(mu.withlink))
  phi <- switch(phi.link, "fishZ" = invfishZ(phi.withlink), "logit" = invlogit(phi.withlink))

  val <- sum(VGAM::dbetabinom(W, M, prob = mu, rho = phi, log = TRUE))

  return(val)
}
