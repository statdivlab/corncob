#' Parameter Gradient Vector
#'
#' Used for internal optimization. Not intended for users.
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
#' @return Gradient of likelihood with respect to parameters
gr_full <- function(theta, W, M, X, X_star, np, npstar, link, phi.link, logpar = TRUE) {
  # extract matrix of betas (np x 1), first np entries
  b      <- utils::head(theta, np)
  # extract matrix of beta stars (npstar x 1), last npstar entries
  b_star <- utils::tail(theta, npstar)

  mu.withlink <- X %*% b
  phi.withlink <- X_star %*% b_star
  mu <- switch(link, "logit" = invlogit(mu.withlink))
  phi <- switch(phi.link, "fishZ" = invfishZ(phi.withlink), "logit" = invlogit(phi.withlink))

  gam <- phi/(1 - phi)

  # Hold digammas
  dg1 <- digamma(M - (mu + W * gam - 1)/gam)
  dg2 <- digamma((1 - mu)/gam)
  dg3 <- digamma(mu/gam)
  dg4 <- digamma(mu/gam + W)
  dg5 <- digamma(1/gam)
  dg6 <- digamma(M + 1/gam)

  # Hold partials - this part is fully generalized
  dldmu <- (-dg1 + dg2 - dg3 + dg4)/gam
  dldgam <- (-dg5 + dg6 + (mu - 1) * (dg1 - dg2) + mu * (dg3 - dg4))/(gam^2)

  # NOTE: Below depends on link functions! Right now for logit and fishZ
  tmp_b <- switch(link, "logit" = mu * (1 - mu) * dldmu)
  tmp_bstar <- switch(phi.link, "fishZ" = (gam + 0.5) * dldgam, "logit" = gam * dldgam)

  # Add in covariates
  g_b <- c(crossprod(tmp_b, X))
  g_bstar <- c(crossprod(tmp_bstar, X_star))

  return(-c(g_b, g_bstar))
}
