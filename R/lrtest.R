#' Likelihood ratio test
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null an object of class \code{bbdml}, should be nested within \code{mod}
#'
#' @return P-value from likelihood ratio test.
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
#' lrtest(mod1, mod2)
#' }
#' @export
lrtest <- function(mod, mod_null) {
  checkNested(mod, mod_null)
  dof.dif <- mod$df.model - mod_null$df.model
  chi.val <- 2 * abs(mod$logL - mod_null$logL)
  pvalue <- stats::pchisq(chi.val, dof.dif, lower.tail = FALSE)
  return(pvalue)
}
