#' Compute sandwich standard error
#'
#' @param mod model fit from bbdml
#' @param numerical Boolean indicator for numerical score and hessian. Not as stable. Defaults to FALSE
#'
#' @return Sandwich variance-covariance matrix
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#'
#' @export
sandSE <- function(mod, numerical = FALSE) {
  # Form A^-1 * B * A^-1
  A <- hessian(mod, numerical = numerical)
  B <- score(mod, numerical = numerical, forHess = TRUE)
  Ainv <- solve(A)
  return(Ainv %*% B %*% Ainv)
}
