#' Compute score at the MLE
#'
#' @param mod an object of class \code{bbdml}
#' @param numerical Boolean. Defaults to \code{FALSE}. Indicator of whether to use the numeric Hessian and score (not recommended).
#' @param forHess Boolean. Defaults to \code{FALSE}. Bryan: Indicator of whether to put in vector form. This parameter is not intended for users. Amy: actually returns robust estimate of variance of score: $hat{B}(hat{theta}) = sum_i G(hat{theta}; W_i) G(hat{theta}; W_i)^T $.
#'
#' @return Score at the MLE. For $G(theta, w)$ score function, returns $sum_i G(hat{theta}, W_i)$ if forHess = FALSE.
#'
#' @examples
#' data(soil_phylum_small)
#' mod <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small)
#' score(mod)
#'
#' @export
score <- function(mod, numerical = FALSE, forHess = FALSE) {
  mu <- mod$mu.resp
  phi <- mod$phi.resp
  W <- mod$W
  M <- mod$M
  X <- mod$X.mu
  X_star <- mod$X.phi
  link <- mod$link
  phi.link <- mod$phi.link

  npx <- ncol(X)
  npw <- ncol(X_star)

  if (numerical) {
    ## care needed if trying to calculate at somewhere other than MLE
    stopifnot(length(mod$param) == ncol(mod$X.mu) + ncol(mod$X.phi))

    return(numDeriv::grad(func = dbetabin, x = mod$param, W = W, M = M,
                          X = X, X_star = X_star, np = npx, npstar = npw,
                          link = mod$link, phi.link = mod$phi.link))
  }

  ### STEP 2 - Gradient

  # define gam
  gam <- phi/(1 - phi)

  # Hold digammas
  dg1 <- digamma(1/gam)
  dg2 <- digamma(M + 1/gam)
  dg3 <- digamma(M - (mu + W * gam - 1)/gam)
  dg4 <- digamma((1 - mu)/gam)
  dg5 <- digamma(mu/gam)
  dg6 <- digamma(mu/gam + W)

  # Hold partials - this part is fully generalized
  dldmu <- (-dg3 + dg4 - dg5 + dg6)/gam
  dldgam <- (-dg1 + dg2 + (mu - 1) * (dg3 - dg4) + mu *
               (dg5 - dg6))/gam^2

  # NOTE: Below depends on link functions! Right now for logit and fishZ
  tmp_b <- switch(link, "logit" = mu * (1 - mu) * dldmu)
  tmp_bstar <- switch(phi.link, "fishZ" = (gam + 0.5) * dldgam, "logit" = gam * dldgam)


  # Add in covariates
  g_b <- c(crossprod(tmp_b, X))
  g_bstar <- c(crossprod(tmp_bstar, X_star))

  # Keep in terms of subject-specific score, useful for sandwich
  if (forHess) {
    V <- matrix(0, nrow = npx + npw, ncol = npx + npw)
    nu <- rep(0, npx + npw)
    # sample size
    n <- length(W)
    for (i in 1:n) {
      x <- X[i,]
      b <- tmp_b[i,]
      x_star <- X_star[i,]
      bst <- tmp_bstar[i,]
      nu <- c(x * b, x_star * bst)
      V <- V + tcrossprod(nu)
    }
    return(V)
  }

  return(c(g_b, g_bstar))
}
