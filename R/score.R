#' Compute score
#'
#' @param mod model fit from bbdml
#' @param numerical Boolean numerical score. Not as stable. Defaults to FALSE
#' @param forHess Boolean for whether to use to approximate Hessian. Defaults to FALSE.
#'
#' @return Analytic score
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#'
#' @export
score <- function(mod, numerical = FALSE, forHess = FALSE) {
  mu <- mod$mu.resp
  phi <- mod$phi.resp
  gam <- phi/(1 - phi)
  N <- mod$W
  M <- mod$M
  X <- mod$X.mu
  W <- mod$X.phi
  link <- mod$link
  phi.link <- mod$phi.link

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
  # n by p
  if (link == "logit") {
    tmp_b <- mu * (1 - mu) * dldmu
  }
  if (phi.link == "fishZ") {
    tmp_bstar <- (gam + 0.5) * dldgam
  } else if (phi.link == "logit") {
    tmp_bstar <- gam * dldgam
  }


  # Keep in terms of subject-specific score, useful for sandwich
  if (forHess) {
    V <- matrix(0, nrow = npx + npw, ncol = npx + npw)
    nu <- rep(0, npx + npw)
    # sample size
    n <- length(N)
    for (i in 1:n) {
      x <- X[i,]
      b <- tmp_b[i,]
      w <- W[i,]
      bst <- tmp_bstar[i,]
      nu <- c(x * b, w * bst)
      V <- V + tcrossprod(nu)
    }
    return(V)
  }

  # Add in covariates
  g_b <- c(crossprod(tmp_b, X))
  g_bstar <- c(crossprod(tmp_bstar, W))

  return(c(g_b, g_bstar))
}
