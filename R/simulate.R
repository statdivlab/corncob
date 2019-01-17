#' Simulate from beta-binomial model
#'
#' @param mod an object of class \code{bbdml}
#' @param nsim Integer. Number of simulations
#' @param seed Optional integer to set a random seed
#' @param ... There are no additional parameters at this time.
#'
#' @return \code{nsim} simulations from \code{mod}
#'
#' @importFrom stats simulate
#' @export
simulate.bbdml <- function(object, nsim, seed = NULL, ...) {
  if (!is.null(seed)) {set.seed(seed)}
  return(VGAM::rbetabinom(n = nsim, size = object$M, prob = object$mu.resp, rho = object$phi.resp))
}
