#' Likelihood ratio test
#'
#' @param mod1 model fit from bbdml, nested within \code{mod2}
#' @param mod2 model fit from bbdml
#'
#' @return Matrix with likelihood ratio test statistic and p-value
#'
#' @examples
#' \dontrun{
#' TODO
#' }
lrtest <- function(mod1, mod2) {
  dof.dif <- mod2$df.model - mod1$df.model
  chi.val <- 2 * abs(mod2$logL - mod1$logL)
  pvalue <- stats::pchisq(chi.val, dof.dif, lower.tail = FALSE)
  return(pvalue)
}
