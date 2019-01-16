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
    warning("Singular Hessian! Cannot calculate p-values in this setting.")
    np <- length(mod$param)
    se <- tvalue <- pvalue <- rep(NA, np)
  } else {
    # Standard errors
    se <- sqrt(diag(covMat))
    # test statistic
    tvalue <- mod$param/se
    # P-value
    pvalue <- 2*stats::pt(-abs(tvalue), mod$df.residual)
  }
  # make table
  coef.table <- cbind(mod$param, se, tvalue, pvalue)
  dimnames(coef.table) <- list(names(mod$param),
                               c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  return(coef.table)
}





#' Wald-type chi-squared test statistic
#'
#' @param mod unrestricted model fit from \code{bbdml}
#' @param restrictions Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the abundance to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the abundance.
#' @param restrictions.phi Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the dispersion to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the dispersion.
#'
#' @return Test statistic.
waldchisq_test <- function(mod, restrictions = NULL, restrictions.phi = NULL) {
  if (length(restrictions) == 0 && length(restrictions.phi) == 0) {
    stop("No restrictions provided!")
  }
  # Covariance matrix - I_n^-1
  covMat <- try(chol2inv(chol(hessian(mod))), silent = TRUE)
  if (class(covMat) == "try-error") {
    stop("Singular Hessian!")
  }

  #### Get index of terms to be tested ####
  if (!is.null(restrictions) && !is.numeric(restrictions)) {
    if (is.character(restrictions)) {
      restrictions <- getRestrictionTerms(mod = mod, restrictions = restrictions)$mu
    } else {
      stop("restrictions must be either character vector or integer vector!")
    }
  }
  if (!is.null(restrictions.phi) && !is.numeric(restrictions.phi)) {
    if (is.character(restrictions.phi)) {
      restrictions.phi <- getRestrictionTerms(mod = mod, restrictions.phi = restrictions.phi)$phi
    } else {
      stop("restrictions must be either character vector or integer vector!")
    }
  }
  if (is.null(attr(restrictions.phi, "added"))) {
    restrictions.phi <- restrictions.phi + mod$np.mu
  }

  index <- c(restrictions, restrictions.phi)


  cov_test <- covMat[index, index]
  par_test <- mod$param[index]
  chi.val <- c(crossprod(par_test, chol2inv(chol(cov_test))) %*% par_test)
  attr(chi.val, "df") <- length(index)
  return(chi.val)
}

#' Wald-type chi-squared test
#'
#' @param mod unrestricted model fit from \code{bbdml}
#' @param mod_null Optional. Restricted model fit from \code{bbdml}. If not included, need to include \code{restrictions} or \code{restrictions.phi}.
#' @param restrictions Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the abundance to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the abundance.
#' @param restrictions.phi Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the dispersion to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the dispersion.
#'
#' @return P-value.
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
waldchisq <- function(mod, mod_null = NULL, restrictions = NULL, restrictions.phi = NULL) {
  if (!is.null(mod_null)) {
    tmp <- getRestrictionTerms(mod = mod, mod_null = mod_null)
    restrictions <- tmp$mu
    restrictions.phi <- tmp$phi
  }
  chi.val <- try(waldchisq_test(mod, restrictions = restrictions, restrictions.phi = restrictions.phi),
                 silent = TRUE)
  if (class(chi.val) == "try-error") {
    return(NA)
  }
  dof.dif <- attr(chi.val, "df")
  return(stats::pchisq(chi.val, dof.dif, lower.tail = FALSE))
}

