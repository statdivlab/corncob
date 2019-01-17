#' Wald-type t test
#'
#' @param mod an object of class \code{bbdml}
#'
#' @return Matrix with wald test statistics and p-values. Only performs univariate tests.
#'
#' @examples
#' \dontrun{
#' data(soil_phylo)
#' soil <- soil_phylo %>%
#' phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
#' tax_glom("Phylum")
#' mod1 <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil)
#' waldt(mod1)
#' }
#' @export
waldt <- function(mod) {
  # Covariance matrix
  covMat <- try(chol2inv(chol(hessian(mod))), silent = TRUE)
  if (class(covMat) == "try-error") {
    warning("Singular Hessian! Cannot calculate p-values in this setting.", .immediate = TRUE)
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
#' This is a helper function and not intended for users
#'
#' @param mod an object of class \code{bbdml}
#' @param restrictions Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the abundance to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the abundance.
#' @param restrictions.phi Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the dispersion to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the dispersion.
#'
#' @return Test statistic for Wald test.
#'
#' @examples
#' \dontrun{
#' data(soil_phylo)
#' soil <- soil_phylo %>%
#' phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
#' tax_glom("Phylum")
#' mod1 <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil)
#'
#' waldchisq_test(mod = mod1, restrictions = 2)
#' }
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
      stop("restrictions.phi must be either character vector or integer vector!")
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
#' @param mod an object of class \code{bbdml}
#' @param mod_null Optional. An object of class \code{bbdml}, should be nested within \code{mod}. If not included, need to include \code{restrictions} or \code{restrictions.phi}.
#' @param restrictions Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the abundance to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the abundance.
#' @param restrictions.phi Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the dispersion to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the dispersion.
#'
#' @return P-value from Wald test.
#'
#' @examples
#' \dontrun{
#' data(soil_phylo)
#' soil <- soil_phylo %>%
#' phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
#' tax_glom("Phylum")
#' mod1 <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil)
#'
#' mod2 <- bbdml(formula = OTU.1 ~ 1,
#' phi.formula = ~ 1,
#' data = soil)
#'
#' # Example using mod_null
#' waldchisq(mod = mod1, mod_null = mod2)
#'
#' # Example using restrictions and restrictions.phi
#' waldchisq(mod = mod1, restrictions = 2, restrictions.phi = 2)
#' waldchisq(mod = mod1, restrictions = "DayAmdmt", restrictions.phi = "DayAmdmt")
#' waldchisq(mod = mod1, restrictions = 2, restrictions.phi = "DayAmdmt")
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

