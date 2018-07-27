#' Wald test
#'
#' @param mod model fit from bbdml
#'
#' @return Matrix with wald test statistics and p-values
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
waldtest <- function(mod) {
  # Covariance matrix
  covMat <- try(chol2inv(chol(hessian(mod))), silent = TRUE)
  if (class(covMat) == "try-error") {
    stop("Singular Hessian!")
  }
  # Standard errors
  se <- sqrt(diag(covMat))
  # test statistic
  tvalue <- mod$param/se
  # P-value
  pvalue <- 2*stats::pt(-abs(tvalue), mod$df.residual)
  # make table
  coef.table <- cbind(mod$param, se, tvalue, pvalue)
  dimnames(coef.table) <- list(names(mod$param),
                               c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  return(coef.table)
}
