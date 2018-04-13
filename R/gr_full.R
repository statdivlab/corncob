#' Parameter Gradient Vector
#'
#' @param theta parameters
#' @param W absolute abundance
#' @param M sample size
#' @param X mean covariates
#' @param X_star overdispersion covariates
#' @param np number of mean parameters
#' @param npstar number of overdisperion parameters
#' @param link Link function for mean, defaults to "logit"
#' @param phi.link Link function for overdispersion, defaults to "fishZ"
#' @param logpar Indicator of log-likelihood, defaults to TRUE
#'
#' @return Gradient of likelihood with respect to parameters
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#'
#' @export
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
  tmp_b <- mu * (1 - mu) * dldmu
  tmp_bstar <- (gam + 0.5) * dldgam

  # Add in covariates
  g_b <- c(crossprod(tmp_b, X))
  g_bstar <- c(crossprod(tmp_bstar, X_star))

  return(-c(g_b, g_bstar))
}
