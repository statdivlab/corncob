#' Likelihood ratio test
#'
#' @param mod model fit from bbdml
#' @param mod_null model fit from bbdml, nested within \code{mod}
#'
#' @return P-value.
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
lrtest <- function(mod, mod_null) {
  checkNested(mod, mod_null)
  dof.dif <- mod$df.model - mod_null$df.model
  chi.val <- 2 * abs(mod$logL - mod_null$logL)
  pvalue <- stats::pchisq(chi.val, dof.dif, lower.tail = FALSE)
  return(pvalue)
}
