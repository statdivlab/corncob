#' Simulate from beta-binomial model
#'
#' @param object an object of class \code{bbdml}
#' @param nsim Integer. Number of simulations
#' @param seed Optional integer to set a random seed
#' @param ... There are no additional parameters at this time.
#'
#' @return \code{nsim} simulations from \code{object}
#'
#' @importFrom stats simulate
#' @export
simulate.bbdml <- function(object, nsim, seed = NULL, ...) {
  if (!is.null(seed)) {set.seed(seed)}
  if (!all(round(object$M) == object$M)) {
    stop(paste("Can't simulate from a beta-binomial distribution without integer M.\n",
               "If this is for bootstrap testing, please choose a different testing approach."))
  }
  return(VGAM::rbetabinom(n = nsim, size = object$M, prob = object$mu.resp, rho = object$phi.resp))
}

