#' Simulate from BBD model fit
#'
#' @param object model output from \code{\link{bbdml}}
#' @param nsim number of simulations
#' @param seed object specifying random seed if desired
#' @param ... optional additional arguments
#'
#' @return nsim simulations from mod
#'
#' @importFrom stats simulate
#' @export
simulate.bbdml <- function(object, nsim, seed = NULL, ...) {
  if (!is.null(seed)) {set.seed(seed)}
  return(VGAM::rbetabinom(n = nsim, size = object$M, prob = object$mu.resp, rho = object$phi.resp))
}
