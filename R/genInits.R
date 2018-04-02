#' Generate initialization for optimization
#'
#' @param nstart Number of starts for optimization
#' @param max.start.time Maximum amount of time to select each start in seconds
#' @param W absolute abundance
#' @param M sample size
#' @param X mean covariates
#' @param X_star overdispersion covariates
#' @param np number of mean parameters
#' @param npstar number of overdisperion parameters
#' @param link Link function for mean
#' @param phi.link Link function for overdispersion
#' @param logpar Indicator of log-likelihood
#' @param lower Lower bound for parameter values. Defaults to -20
#' @param upper Upper bound for parameter values. Defaults to 20.
#'
#' @return Matrix of initializations
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
genInits <- function(nstart, max.start.time,
                     W, M,
                     X, X_star,
                     np, npstar,
                     link, phi.link, logpar,
                     lower = NULL, upper = NULL) {

  ll_init <- function(theta, W, M, X, X_star, np, npstar, link, phi.link) {
    # extract matrix of betas (np x 1), first np entries
    b      <- utils::head(theta, np)
    # extract matrix of beta stars (npstar x 1), last npstar entries
    b_star <- utils::tail(theta, npstar)

    mu.withlink <- X %*% b
    phi.withlink <- X_star %*% b_star
    mu <- switch(link, "logit" = invlogit(mu.withlink))
    phi <- switch(phi.link, "fishZ" = invfishZ(phi.withlink))

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
    return(value)
  }

  if (is.null(lower)) {
    lower <- rep(-20, np + npstar)
  }
  if (is.null(upper)) {
    upper <- rep(20, np + npstar)
  }

  inits <- matrix(NA, nrow = nstart, ncol = np + npstar)
  for (i in 1:nstart) {
    inits[i,] <- GenSA::GenSA(fn = ll_init, lower = lower, upper = upper, W = W, M = M, X = X, X_star = X_star, np = np,
                         npstar = npstar, link = link, phi.link = phi.link, control = list(max.time = max.start.time))$par
  }


  return(inits)
}
