#' Get index of restricted terms for Wald test
#'
#' Created as a convenient helper function. Not intended for users.
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null Optional. An object of class \code{bbdml}. Defaults to \code{NULL}
#' @param restrictions Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the abundance to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the abundance.
#' @param restrictions.phi Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the dispersion to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the dispersion.
#'
#' @return A list with \code{mu} representing the index of the restricted covariates associated with abundance and \code{phi} representing the index of the restricted covarates associated with the dispersion
getRestrictionTerms <- function(mod, mod_null = NULL, restrictions = NULL, restrictions.phi = NULL) {

  # Term labels
  full.mu <- attr(stats::terms(mod$formula), "term.labels")
  full.phi <- attr(stats::terms(mod$phi.formula), "term.labels")

  # Term assignments
  assigns.mu <- attr(mod$X.mu, "assign")
  assigns.phi <- attr(mod$X.phi, "assign")

  sortInteraction <- function(x) {
    # Sorts interaction effects so they always match
    sapply(lapply(strsplit(x, ":"), sort), paste, collapse = ":")
  }

  if (!is.null(mod_null)) {
    checkNested(mod, mod_null)

    restrict.mu <- attr(stats::terms(mod_null$formula), "term.labels")
    restrict.phi <- attr(stats::terms(mod_null$phi.formula), "term.labels")

    restrictions <- setdiff(full.mu, restrict.mu)
    restrictions.phi <- setdiff(full.phi, restrict.phi)

    # Get index combining
    index.mu <- which(assigns.mu %in% match(sortInteraction(restrictions),
                                            sortInteraction(full.mu)))
    index.phi <- which(assigns.phi %in% match(sortInteraction(restrictions.phi),
                                              sortInteraction(full.phi)))

    # Add intercept term if applicable
    if (attr(stats::terms(mod$formula), "intercept") == 1 &&
        attr(stats::terms(mod_null$formula), "intercept") == 0) {
      index.mu <- c(1, index.mu)
    }

    if (attr(stats::terms(mod$phi.formula), "intercept") == 1 &&
        attr(stats::terms(mod_null$phi.formula), "intercept") == 0) {
      index.phi <- c(1, index.phi)
    }

    index.phi <- index.phi + mod$np.mu
    attr(index.phi, "added") <- TRUE

    if (length(index.mu) == 0) {
      index.mu <- NULL
    }

    if (length(index.phi) == 0) {
      index.phi <- NULL
    }

    return(list(mu = index.mu, phi = index.phi))
  }


  if (is.numeric(restrictions) || is.null(restrictions)) {
    index.mu <- restrictions
  }

  if (is.null(restrictions.phi)) {
    index.phi <- NULL
  }


  if (length(restrictions) != 0 && !is.numeric(restrictions)) {
    if (is.character(restrictions)) {
      index.mu <- which(assigns.mu %in% match(sortInteraction(restrictions),
                                              sortInteraction(full.mu)))
      if ("(Intercept)" %in% restrictions) {
        index.mu <- c(1, index.mu)
      }
    } else {
      stop("restrictions must be either character vector or integer vector!")
    }
  }


  if (is.numeric(restrictions.phi) && is.null(attr(restrictions.phi, "added"))) {
    index.phi <- restrictions.phi + mod$np.mu
    attr(index.phi, "added") <- TRUE
  }

  if (length(restrictions.phi) != 0 && !is.numeric(restrictions.phi)) {
    if (is.character(restrictions.phi)) {
      index.phi <- which(assigns.phi %in% match(sortInteraction(restrictions.phi),
                                                sortInteraction(full.phi)))
      if ("(Intercept)" %in% restrictions.phi) {
        index.phi <- c(1, index.phi)
      }
      index.phi <- index.phi + mod$np.mu
      attr(index.phi, "added") <- TRUE
    } else {
      stop("restrictions.phi must be either character vector or integer vector!")
    }
  }

  if (length(index.mu) == 0) {
    index.mu <- NULL
  }

  if (length(index.phi) == 0) {
    index.phi <- NULL
  }

  return(list(mu = index.mu, phi = index.phi))
}
