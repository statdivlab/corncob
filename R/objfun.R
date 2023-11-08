#' Objective function
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
#'
#' @return List of negative log-likelihood, gradient, and hessian
objfun <- function(theta, W, M, X, X_star, np, npstar, link, phi.link) {

  ### STEP 1 - Negative Log-likelihood

  # extract matrix of betas (np x 1), first np entries
  b      <- utils::head(theta, np)
  # extract matrix of beta stars (npstar x 1), last npstar entries
  b_star <- utils::tail(theta, npstar)

  mu.withlink <- X %*% b
  phi.withlink <- X_star %*% b_star
  mu <- switch(link, "logit" = invlogit(mu.withlink))
  phi <- switch(phi.link, "fishZ" = invfishZ(phi.withlink), "logit" = invlogit(phi.withlink))

  val <- suppressWarnings(sum(dbetabinom_cts(W, M, prob = mu, rho = phi, log = TRUE)))
  if (is.nan(val) || any(phi <= sqrt(.Machine$double.eps)) || any(phi >= 1 - sqrt(.Machine$double.eps))) {
    return(list(value = Inf))
  }
  value <- -val

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

  gradient <- -c(g_b, g_bstar)

  ### STEP 3 - Hessian

  tg1 <- trigamma(M - (mu + W * gam - 1)/gam)
  tg2 <- trigamma((1 - mu)/gam)
  tg3 <- trigamma(mu/gam)
  tg4 <- trigamma(mu/gam + W)
  tg5 <- trigamma(1/gam)
  tg6 <- trigamma(M + 1/gam)

  dldmu2 <- (tg1 - tg2 - tg3 + tg4)/gam^2
  dldgam2 <- (2 * gam * dg1 + tg5 - 2 * gam * dg2 - tg6 +
                (mu - 1)^2 * tg1 - 2 * gam * (mu - 1) * dg3 - mu^2 *
                tg3 + mu^2 * tg4 + (mu - 1)^2 * (-tg2) + 2 * gam * (mu -
                                                                  1) * dg4 - 2 * gam * mu * dg5 + 2 * gam * mu * dg6)/gam^4
  dldmdg <- (gam * (dg3 - dg4 + dg5 - dg6) + (mu - 1) * (tg2 -
                                                        tg1) + mu * (tg3 - tg4))/gam^3

  dpdb <- switch(link, logit = t(X * c(mu * (1 - mu))))
  dgdb <- switch(phi.link, fishZ = t(X_star * c(gam + 0.5)),
                 logit = t(X_star * c(gam)))

  mid4 <- switch(link, logit = c(dldmu * mu * (1 - mu) * (1 - 2 * mu)))
  mid5 <- switch(phi.link, fishZ = c(dldgam * (gam + 0.5)),
                 logit = c(dldgam * gam))
  term4 <- crossprod(X, diag(mid4)) %*% X
  term5 <- crossprod(X_star, diag(mid5)) %*% X_star
  term1 <- dpdb %*% tcrossprod(diag(c(-dldmu2)), dpdb)
  term2 <- dpdb %*% tcrossprod(diag(c(-dldmdg)), dgdb)
  term3 <- dgdb %*% tcrossprod(diag(c(-dldgam2)), dgdb)

  hessian <- cbind(rbind(term1 - term4, t(term2)),rbind(term2, term3 - term5))

  return(list(value = value, gradient = gradient, hessian = hessian))

}
