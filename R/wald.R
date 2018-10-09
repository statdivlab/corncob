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





#' Wald-type chi-squared test statistic
#'
#' @param mod unrestricted model fit from \code{bbdml}
#' @param restrictions Numeric vector indicating the parameters to test, or character vector with name of variable to test.
#' @param testonly Optional. If set to \code{"mu"}, then only the effect of the variable provided in \code{restrictions} on the mean is tested. If set to \code{"phi"}, then only the effect of the variable provided in \code{restrictions} on the variance is tested. Otherwise, the effects are jointly tested.
#'
#' @return Test statistic.
waldchisq_test <- function(mod, restrictions, testonly = NULL) {
  # Covariance matrix - I_n^-1
  covMat <- try(chol2inv(chol(hessian(mod))), silent = TRUE)
  if (class(covMat) == "try-error") {
    stop("Singular Hessian!")
  }

  if (length(restrictions) == 0) {
    stop("restrictions is of length 0!")
  }
  #### Get index of terms to be tested ####
  if (is.character(restrictions)) {
    # Term labels
    allterms.mu <- attr(stats::terms(stats::model.frame(mod$formula, data = mod$dat)), "term.labels")
    allterms.phi <- attr(stats::terms(stats::model.frame(mod$phi.formula, data = mod$dat)), "term.labels")
    # Term assignments
    assigns.mu <- attr(mod$X.mu, "assign")
    assigns.phi <- attr(mod$X.phi, "assign")
    sortInteraction <- function(x){
      # Sorts interaction effects so they always match
      sapply(lapply(strsplit(x,":"), sort), paste, collapse = ":")
    }
    # Get index combining
    index.mu <- which(assigns.mu %in% match(sortInteraction(restrictions), sortInteraction(allterms.mu)))
    index.phi <- which(assigns.phi %in% match(sortInteraction(restrictions), sortInteraction(allterms.phi))) + mod$np.mu
    if (is.null(testonly)) {
      index <- c(index.mu, index.phi)
    } else {
      if (testonly == "mu") {
        index <- index.mu
      } else if (testonly == "phi") {
        index <- index.phi
      }
    }
  } else if (is.numeric(restrictions)) {
    index <- restrictions
  } else {
    stop("restrictions must be either character vector or integer vector!")
  }

  cov_test <- covMat[index, index]
  par_test <- mod$param[index]
  chi.val <- crossprod(par_test, chol2inv(chol(cov_test))) %*% par_test
  attr(chi.val, "df") <- length(index)
  return(chi.val)
}

#' Wald-type chi-squared test
#'
#' @param mod unrestricted model fit from \code{bbdml}
#' @param mod_null Optional. Restricted model fit from \code{bbdml}. If not included, need to include \code{restrictions}.
#' @param restrictions Optional. Numeric vector indicating the parameters to test, or character vector with name of variable to test.
#' @param testonly Optional. Defaults to \code{NULL}. If set to \code{"mu"}, then only the effect of the variable provided in \code{restrictions} on the mean is tested. If set to \code{"phi"}, then only the effect of the variable provided in \code{restrictions} on the variance is tested. Otherwise, the effects are jointly tested.
#'
#' @return P-value.
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
waldchisq <- function(mod, mod_null = NULL, restrictions = NULL, testonly = NULL) {
  if (is.null(restrictions)) {
    restrictions <- getRestrictionTerms(mod = mod, mod_null = mod_null)
    testonly <- attr(restrictions, "testonly")
  }
  chi.val <- waldchisq_test(mod, restrictions, testonly)
  dof.dif <- attr(chi.val, "df")
  return(stats::pchisq(chi.val, dof.dif, lower.tail = FALSE))
}

