#' Get index of restricted terms for Wald test
#'
#' @param mod unrestricted model fit from \code{bbdml}
#' @param mod_null restricted model fit from \code{bbdml}
#'
#' @return Numeric vector of index of restricted terms.
getRestrictionTerms <- function(mod, mod_null) {
  altterms.mu <- attr(stats::terms(stats::model.frame(mod$formula, data = mod$dat)), "term.labels")
  altterms.phi <- attr(stats::terms(stats::model.frame(mod$phi.formula, data = mod$dat)), "term.labels")
  nullterms.mu <- attr(stats::terms(stats::model.frame(mod_null$formula, data = mod_null$dat)), "term.labels")
  nullterms.phi <- attr(stats::terms(stats::model.frame(mod_null$phi.formula, data = mod_null$dat)), "term.labels")

  res_mu <- setdiff(altterms.mu, nullterms.mu)
  res_phi <- setdiff(altterms.phi, nullterms.phi)
  testonly <- NULL
  if (length(res_mu) == 0) {
    testonly <- "phi"
  } else if (length(res_phi) == 0) {
    testonly <- "mu"
  }
  restrictions <- union(res_mu, res_phi)
  attr(restrictions, "testonly") <- testonly
  return(restrictions)
}
