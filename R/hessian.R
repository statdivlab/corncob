#' Compute Hessian
#'
#' @param mod model fit from bbdml
#' @param numerical Boolean numerical Hessian. Not as stable. Defaults to FALSE
#'
#' @return Analytic Hessian
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#'
#' @export
hessian <- function(mod, numerical = FALSE) {
  mu <- mod$mu.resp
  phi <- mod$phi.resp
  gam <- phi/(1 - phi)
  N <- mod$W
  M <- mod$M
  X <- mod$X.mu
  W <- mod$X.phi


  npx <- ncol(X)
  npw <- ncol(W)
  H <- matrix(0, nrow = npx + npw, ncol = npx + npw)

  if (numerical) {
    return(numDeriv::hessian(func = dbetabin_pos, x = mod$param, W = N, M = M,
                             X = X, X_star = W, np = npx, npstar = npw,
                             link = mod$link, phi.link = mod$phi.link))
  }

  # This part is fully generalized
  for (i in 1:length(mu)) {
    m <- mu[i]
    g <- gam[i]
    y <- N[i]
    n <- M[i]
    x <- X[i,]
    w <- W[i,]
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

    dpdb2 <- as.matrix(Matrix::bdiag(dpdb2, matrix(0, nrow = npw, ncol = npw)))
    dgdb2 <- as.matrix(Matrix::bdiag(matrix(0, nrow = npx, ncol = npx), dgdb2))

    term1 <- (-dldmu2) * tcrossprod(dpdb)
    term2 <- (-dldmdg) * (tcrossprod(dpdb, dgdb) + tcrossprod(dgdb, dpdb))
    term3 <- (-dldgam2) * tcrossprod(dgdb)
    term4 <- dldmu*dpdb2
    term5 <- dldgam*dgdb2
    H <- H + term1 + term2 + term3 - term4 - term5
  }
  return(H)
}
