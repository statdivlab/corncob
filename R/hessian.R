#' Compute Hessian matrix at the MLE
#'
#' @param mod an object of class \code{bbdml}
#' @param numerical Boolean. Defaults to \code{FALSE}. Indicator of whether to use the numeric Hessian (not recommended).
#' @return Hessian matrix at the MLE
#'
#' @examples
#' data(soil_phylum_small)
#' mod <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small)
#' hessian(mod)
#'
#' @export
hessian <- function(mod, numerical = FALSE) {
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
    stopifnot(length(mod$param) != ncol(mod$X.mu) + ncol(mod$X.phi))

    return(numDeriv::hessian(func = dbetabin_neg, x = mod$param, W = W, M = M,
                             X = X, X_star = X_star, np = npx, npstar = npw,
                             link = mod$link, phi.link = mod$phi.link))
  }

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

  H <- cbind(rbind(term1 - term4, t(term2)),rbind(term2, term3 - term5))
  return(H)
}
