#' Get index of restricted terms for Wald test
#'
#' @param mod unrestricted model fit from \code{bbdml}
#' @param mod_null restricted model fit from \code{bbdml}. Defaults to \code{NULL}.
#' @param restrictions Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the abundance to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the abundance.
#' @param restrictions.phi Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the dispersion to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the dispersion.
#'
#' @return A numeric vector of index of restricted terms.
getRestrictionTerms <- function(mod, mod_null = NULL, restrictions = NULL, restrictions.phi = NULL) {
  if (!is.null(mod_null)) {
    checkNested(mod, mod_null)
    restrictions <- substring(setdiff(names(mod$b.mu), names(mod_null$b.mu)), 4)
    restrictions.phi <- substring(setdiff(names(mod$b.phi), names(mod_null$b.phi)), 5)
  }
  if (length(restrictions) != 0 && !is.numeric(restrictions)) {
    if (is.character(restrictions)) {
      restrictions <- paste("mu", restrictions, sep = ".")
      restrictions <- which(!is.na(match(names(mod$param), restrictions)))
    } else {
      stop("restrictions must be either character vector or integer vector!")
    }
  }

  if (is.numeric(restrictions.phi) && attr(restrictions.phi, "added") != TRUE) {
    restrictions.phi <- restrictions.phi + mod$np.mu
    attr(restrictions.phi, "added") <- TRUE
  }

  if (length(restrictions.phi) != 0 && !is.numeric(restrictions.phi)) {
    if (is.character(restrictions.phi)) {
      restrictions.phi <- paste("phi", restrictions.phi, sep = ".")
      restrictions.phi <- which(!is.na(match(names(mod$param), restrictions.phi)))
      attr(restrictions.phi, "added") <- TRUE
    } else {
      stop("restrictions must be either character vector or integer vector!")
    }
  }

  if (length(restrictions) == 0) {
    restrictions <- NULL
  }

  if (length(restrictions.phi) == 0) {
    restrictions.phi <- NULL
  }

  return(list(mu = restrictions, phi = restrictions.phi))
}
