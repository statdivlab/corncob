#' Objective function
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

  if (any(mu < 0) || any(mu > 1)) {
    return(list(value = Inf))
  } else if (any(phi <= 0)) {
    return(list(value = Inf))
  } else {
    gam <- phi/(1 - phi)
    a1 <- mu/gam
    a2 <- (1 - mu)/gam
    val <- sum(lbeta(a1 + W, a2 + M - W) - lbeta(a1, a2) + lchoose(M, W))
  }
  value <- -val

  ### STEP 2 - Gradient

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
  if (link == "logit") {
    tmp_b <- mu * (1 - mu) * dldmu
  }
  if (phi.link == "fishZ") {
    tmp_bstar <- (gam + 0.5) * dldgam
  } else if (phi.link == "logit") {
    tmp_bstar <- gam * dldgam
  }



  # Add in covariates
  g_b <- c(crossprod(tmp_b, X))
  g_bstar <- c(crossprod(tmp_bstar, X_star))

  gradient <- -c(g_b, g_bstar)

  ### STEP 3 - Hessian

  H <- matrix(0, nrow = np + npstar, ncol = np + npstar)

  # This part is fully generalized
  for (i in 1:length(mu)) {
    m <- mu[i]
    g <- gam[i]
    y <- W[i]
    n <- M[i]
    x <- X[i,]
    w <- X_star[i,]
    tg1 <- trigamma(n - (m + y*g - 1)/g)
    tg2 <- trigamma((1 - m)/g)
    tg3 <- trigamma(m/g)
    tg4 <- trigamma(m/g + y)

    tg5 <- trigamma(1/g)
    tg6 <- trigamma(n + 1/g)

    dg1 <- digamma(1/g)
    dg2 <- digamma(n + 1/g)
    dg3 <- digamma(n - (m + y*g - 1)/g)
    dg4 <- digamma((1 - m)/g)
    dg5 <- digamma(m/g)
    dg6 <- digamma(m/g + y)

    # Generalizable single derivatives dL
    dldmu <- (-dg3 + dg4 - dg5 + dg6)/g
    dldgam <- (-dg1 + dg2 + (m - 1)*(dg3 - dg4) + m*(dg5 - dg6))/g^2
    # Generalizable double derivatives
    dldmu2 <- (tg1 - tg2 - tg3 + tg4)/g^2
    dldgam2 <- (2*g*dg1 + tg5 - 2*g*dg2 - tg6 + (m - 1)^2*tg1 -
                  2*g*(m - 1)*dg3 - m^2*tg3 + m^2*tg4 + (m - 1)^2*(-tg2) +
                  2*g*(m - 1)*dg4 - 2*g*m*dg5 + 2*g*m*dg6)/g^4
    dldmdg <- (g*(dg3 - dg4 + dg5 - dg6) + (m - 1)*(tg2 - tg1) + m*(tg3 - tg4))/g^3

    # Not generalizeable single dm and dg
    if (link == "logit") {
      dpdb <- x * m * (1 - m)
    }
    if (phi.link == "fishZ") {
      dgdb <- w * (g + 0.5)
    } else if (phi.link == "logit") {
      dgdb <- w * g
    }



    dpdb <- c(dpdb, rep(0, npstar))
    dgdb <- c(rep(0, np), dgdb)
    # Not generalizable double
    if (link == "logit") {
      dpdb2 <- tcrossprod(x) * m * (1 - m) * (1 - 2 * m)
    }
    if (phi.link == "fishZ") {
      dgdb2 <- tcrossprod(w) * (g + 0.5)
    } else if (phi.link == "logit") {
      dgdb2 <- tcrossprod(w) * g
    }



    dpdb2 <- as.matrix(Matrix::bdiag(dpdb2, matrix(0, nrow = npstar, ncol = npstar)))
    dgdb2 <- as.matrix(Matrix::bdiag(matrix(0, nrow = np, ncol = np), dgdb2))

    term1 <- (-dldmu2) * tcrossprod(dpdb)
    term2 <- (-dldmdg) * (tcrossprod(dpdb, dgdb) + tcrossprod(dgdb, dpdb))
    term3 <- (-dldgam2) * tcrossprod(dgdb)
    term4 <- dldmu*dpdb2
    term5 <- dldgam*dgdb2
    H <- H + term1 + term2 + term3 - term4 - term5
  }

  hessian <- H

  return(list(value = value, gradient = gradient, hessian = hessian))




}
