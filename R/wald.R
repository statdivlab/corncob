#' Wald-type t test
#'
#' @param mod model fit from \code{bbdml}
#'
#' @return Matrix with wald test statistics and p-values. Only performs univariate tests
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
waldt <- function(mod) {
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


#' Wald-type chi-squared test
#'
#' @param mod unrestricted model fit from \code{bbdml}
#' @param restrictions Numeric vector indicating the parameters to test
#'
#' @return Matrix with wald test statistics and p-values
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
waldchisq <- function(mod, restrictions) {
  # Covariance matrix - I_n^-1
  covMat <- try(chol2inv(chol(hessian(mod))), silent = TRUE)
  if (class(covMat) == "try-error") {
    stop("Singular Hessian!")
  }
  dof.dif <- length(restrictions)
  cov_test <- covMat[restrictions, restrictions]
  par_test <- mod$param[restrictions]
  chi.val <- crossprod(par_test, chol2inv(chol(cov_test))) %*% par_test
  return(stats::pchisq(chi.val, dof.dif, lower.tail = FALSE))

}
