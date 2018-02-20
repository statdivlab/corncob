#' Compute score
#'
#' @param mod model fit from bbdml
#' @param numerical Boolean numerical score. Not as stable. Defaults to FALSE
#'
#' @return Analytic score
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#'
#' @export
score <- function(mod, numerical = FALSE) {
  mu <- mod$mu.resp
  phi <- mod$phi.resp
  gam <- phi/(1 - phi)
  N <- mod$dat$Wi
  M <- mod$dat$Ni
  X <- mod$X.mu
  W <- mod$X.phi

  npx <- ncol(X)
  npw <- ncol(W)

  if (numerical) {
    return(numDeriv::grad(func = dbetabin_pos, x = mod$param, W = N, M = M,
                          X = X, X_star = W, np = npx, npstar = npw,
                          link = mod$link, phi.link = mod$phi.link))
  }

  # Hold digammas
  dg1 <- digamma(M - (mu + N * gam - 1)/gam)
  dg2 <- digamma((1 - mu)/gam)
  dg3 <- digamma(mu/gam)
  dg4 <- digamma(mu/gam + N)
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

  return(c(g_b, g_bstar))
}
