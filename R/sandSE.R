#' Compute sandwich estimate of variance-covariance matrix
#'
#' @param mod an object of class \code{bbdml}
#' @param numerical Boolean. Defaults to \code{FALSE}. Indicator of whether to use the numeric Hessian and score (not recommended).
#' @return Sandwich variance-covariance matrix. $hat{A}^{-1} hat{B} hat{A}^{-1}.$
#'
#' @examples
#' data(soil_phylum_small)
#' mod <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small)
#' sand_vcov(mod)
#'
#' @export
sand_vcov <- function(mod, numerical = FALSE) {
  # Form A^-1 * B * A^-1
  A_hat <- hessian(mod, numerical = numerical)
  B_hat <- score(mod, numerical = numerical, get_score_covariance = TRUE)
  A_hat_inv <- solve(A_hat)
  return(A_hat_inv %*% B_hat %*% A_hat_inv)
}

#' Compute sandwich standard errors. Legacy function. Use sand_vcov instead.
#'
#' @param mod an object of class \code{bbdml}
#' @param numerical Boolean. Defaults to \code{FALSE}. Indicator of whether to use the numeric Hessian and score (not recommended).
#' @return Sandwich variance-covariance matrix
#'
#' @examples
#' data(soil_phylum_small)
#' mod <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small)
#' sandSE(mod)
#'
#' @export
sandSE <- function(mod, numerical = FALSE) {
  warning("You called sandSE(). Please use sand_vcov() instead. They are the same function but sand_vcov() has more transparent naming.")
  return(sand_vcov(mod, numerical))
}
