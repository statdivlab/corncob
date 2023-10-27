#' Rao-type t test (model-based or robust)
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null an object of class \code{bbdml}, should be nested within \code{mod}
#'
#' @return P-value from likelihood ratio test.
#'
#'
#' @examples
#' data(soil_phylum_small)
#' mod1 <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small)
#'
#' mod2 <- bbdml(formula = OTU.1 ~ 1,
#' phi.formula = ~ 1,
#' data = soil_phylum_small)
#' raotest(mod1, mod2)
#' @export
raotest <- function(mod, mod_null) {

  dof.dif <- mod$df.model - mod_null$df.model

  stopifnot(checkNested(mod, mod_null))
  stopifnot(mod$robust == mod_null$robust)
  stopifnot(dof.dif > 0)

  robust <- mod$robust

  if (robust) {

    stop("Amy hasn't implemented this yet")

    covMat <- sand_vcov(mod)

  } else {

    stop("Amy has implemented this incorrectly!!!")

    #### currently this just calculates the likelihood under the null, at the null
    #### we actually need to calculate the likelihood under the alternative, at the null


    inv_fish_info_null <- try(chol2inv(chol(hessian(mod_null))), silent = TRUE)
    chisq.val <- as.numeric(score(mod_null) %*% (-inv_fish_info_null) %*% score(mod_null))

  }

  if ("try-error" %in% class(inv_fish_info_null)) {
    warning("Singular Hessian! Cannot calculate p-values in this setting.", immediate. = TRUE)
    pvalue <- NA
  } else {
    pvalue <- stats::pchisq(chisq.val, dof.dif, lower.tail = FALSE)
  }

  return(pvalue)
}


#' rao-type chi-squared test statistic (model-based or robust)
#'
#' This is a helper function and not intended for users
#'
#' @param mod an object of class \code{bbdml}
#' @param restrictions Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the abundance to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the abundance.
#' @param restrictions.phi Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the dispersion to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the dispersion.
#' @param contrasts_DA List. Optional. Constructs a contrast matrix. List elements should be characters specifying contrasts in the parameters within \code{formula}. Note that this is only available with \code{"rao"} value for \code{test}.
#' @param contrasts_DV List. Optional. Constructs a contrast matrix. List elements should be characters specifying contrasts in the parameters within \code{phi.formula}. Note that this is only available with \code{"rao"} value for \code{test}.
#' @param robust Should robust standard errors be used? If not, model-based standard arras are used. Logical, defaults to \code{FALSE}.
#'
#' @return Test statistic for rao test.
raochisq_test <- function(mod, restrictions = NULL, restrictions.phi = NULL,
                           contrasts_DA = NULL, contrasts_DV = NULL, robust = FALSE) {
  if (length(restrictions) == 0 && length(restrictions.phi) == 0 &&
      is.null(contrasts_DA) && is.null(contrasts_DV)) {
    stop("No restrictions or contrasts provided!")
  }


  stop("Amy hasn't implemented this yet")

  if (robust) {
    covMat <- sand_vcov(mod)
  } else {

    if (mod$has_noninteger) {
      warning("Your data has non-integer W or M, and you aren't using robust testing. We will let you do this, but you should consider robust testing instead because your data is cannot be assumed to be drawn from a a beta-binomial distribution.")
    }

    # Covariance matrix - I_n^-1
    covMat <- try(chol2inv(chol(hessian(mod))), silent = TRUE)
  }
  if ("try-error" %in% class(covMat)) {
    stop("Singular Hessian!")
  }
  if (is.null(contrasts_DA) && is.null(contrasts_DV)) {

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
    #end restrictions if, begin contrasts
  } else {

    if (!is.null(contrasts_DA)) {
      requireNamespace("limma")
      contr_vec <- suppressWarnings(limma::makeContrasts(contrasts = contrasts_DA, levels = colnames(mod$X.mu)))
      contr_vec <- c(contr_vec, rep(0, mod$np.phi))
    } else if (!is.null(contrasts_DV)) {
      contr_vec <- suppressWarnings(limma::makeContrasts(contrasts = contrasts_DV, levels = colnames(mod$X.phi)))
      contr_vec <- c(rep(0, mod$np.mu), contr_vec)
    }

    par_contr <- crossprod(contr_vec, mod$param)
    cov_contr <- c(crossprod(contr_vec, covMat)) %*% contr_vec
    chi.val <- c(crossprod(par_contr, chol2inv(chol(cov_contr))) %*% par_contr)
    attr(chi.val, "df") <- 1
    return(chi.val)
  }
}


#' rao-type chi-squared test
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null Optional. An object of class \code{bbdml}, should be nested within \code{mod}. If not included, need to include \code{restrictions} or \code{restrictions.phi}.
#' @param restrictions Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the abundance to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the abundance.
#' @param restrictions.phi Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the dispersion to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the dispersion.
#' @param contrasts_DA List. Optional. Constructs a contrast matrix. List elements should be characters specifying contrasts in the parameters within \code{formula}. Note that this is only available with \code{"rao"} value for \code{test}.
#' @param contrasts_DV List. Optional. Constructs a contrast matrix. List elements should be characters specifying contrasts in the parameters within \code{phi.formula}. Note that this is only available with \code{"rao"} value for \code{test}.
#' @param robust Should robust standard errors be used? If not, model-based standard arras are used. Logical, defaults to \code{FALSE}.
#' @return Matrix with rao test statistics and p-values. Only performs univariate tests.
#'
#' @return P-value from rao test.
#'
#' @examples
#' data(soil_phylum_small)
#' mod1 <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small)
#'
#' mod2 <- bbdml(formula = OTU.1 ~ 1,
#' phi.formula = ~ 1,
#' data = soil_phylum_small)
#'
#' # Example using mod_null
#' raochisq(mod = mod1, mod_null = mod2)
#' raochisq(mod = mod1, mod_null = mod2, robust = TRUE)
#'
#' # Example using restrictions and restrictions.phi
#' raochisq(mod = mod1, restrictions = 2, restrictions.phi = 2)
#' raochisq(mod = mod1, restrictions = "DayAmdmt", restrictions.phi = "DayAmdmt")
#' raochisq(mod = mod1, restrictions = 2, restrictions.phi = "DayAmdmt")
#' raochisq(mod = mod1, restrictions = 2, restrictions.phi = 2, robust = TRUE)
#' @export
raochisq <- function(mod, mod_null = NULL, restrictions = NULL,
                      restrictions.phi = NULL,
                      contrasts_DA = NULL, contrasts_DV = NULL,
                      robust = FALSE) {


  stop("Amy hasn't implemented this yet")

  if (!is.null(mod_null)) {
    tmp <- getRestrictionTerms(mod = mod, mod_null = mod_null)
    restrictions <- tmp$mu
    restrictions.phi <- tmp$phi
  }
  chi.val <- try(raochisq_test(mod, restrictions = restrictions,
                                restrictions.phi = restrictions.phi,
                                contrasts_DA = contrasts_DA,
                                contrasts_DV = contrasts_DV,
                                robust = robust),
                 silent = TRUE)
  if (inherits(chi.val, "try-error")) {
    return(NA)
  }
  dof.dif <- attr(chi.val, "df")
  return(stats::pchisq(chi.val, dof.dif, lower.tail = FALSE))
}

