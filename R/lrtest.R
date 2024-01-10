#' Likelihood ratio test
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null an object of class \code{bbdml}, should be nested within \code{mod}
#'
#' @return P-value from likelihood ratio test.
#'
#' @examples
#' data(soil_phylum_small_otu1)
#' mod1 <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small_otu1)
#'
#' mod2 <- bbdml(formula = cbind(W, M - W) ~ 1,
#' phi.formula = ~ 1,
#' data = soil_phylum_small_otu1)
#' lrtest(mod1, mod2)
#' @export
lrtest <- function(mod, mod_null) {
  checkNested(mod, mod_null)
  dof.dif <- mod$df.model - mod_null$df.model
  chi.val <- 2 * abs(mod$logL - mod_null$logL)
  pvalue <- stats::pchisq(chi.val, dof.dif, lower.tail = FALSE)
  return(pvalue)
}
