#' Parametric bootstrap Wald test
#'
#' @param mod Unrestricted model fit from \code{bbdml}.
#' @param mod_null Restricted model fit from \code{bbdml}.
#' @param B Integer. Defaults to \code{1000}. Number of bootstrap iterations.
#'
#' @return P-value.
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
pbWald <- function(mod, mod_null, B = 1000) {
  tmp <- getRestrictionTerms(mod = mod, mod_null = mod_null)
  restrictions <- tmp$mu
  restrictions.phi <- tmp$phi
  t.observed <- try(waldchisq_test(mod, restrictions = restrictions, restrictions.phi = restrictions.phi), silent = TRUE)
  if (class(t.observed) == "try-error") {
    return(NA)
  }

  BOOT <- rep(NA, B)
  for (j in 1:B) {
    #print(j)
    BOOT[j] <- doBoot(mod = mod, mod_null = mod_null, test = "Wald")
  }
  perc.rank <- function(x, y) (1 + sum(stats::na.omit(y) >= x)) / (length(stats::na.omit(y)) + 1)
  p.val <- perc.rank(t.observed, BOOT)
  return(p.val)
}
