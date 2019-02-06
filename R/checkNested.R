#' Check for nested models
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null an object of class \code{bbdml}
#'
#' @return \code{TRUE} if \code{mod_null} is nested within \code{mod}, otherwise it throws an error.
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
#' checkNested(mod1, mod2)
#' }
#' @export
checkNested <- function(mod, mod_null) {

  # Term labels
  full.mu <- attr(stats::terms(mod$formula), "term.labels")
  full.phi <- attr(stats::terms(mod$phi.formula), "term.labels")
  restrict.mu <- attr(stats::terms(mod_null$formula), "term.labels")
  restrict.phi <- attr(stats::terms(mod_null$phi.formula), "term.labels")

  int.mu <- (attr(stats::terms(mod$formula), "intercept") -
               attr(stats::terms(mod_null$formula), "intercept")) >= 0
  int.phi <- (attr(stats::terms(mod$phi.formula), "intercept") -
                 attr(stats::terms(mod_null$phi.formula), "intercept")) >= 0

  if (!all(c(int.mu, all(restrict.mu %in% full.mu)))) {
    stop("Models for abundance are not nested!")
  }

  if (!all(c(int.phi, all(restrict.phi %in% full.phi)))) {
    stop("Models for dispersion are not nested!")
  }

  if (!identical(mod$dat, mod_null$dat)) {
    stop("Models are not fit to the same data!")
  }

  return(TRUE)
}
