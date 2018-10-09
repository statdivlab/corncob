#' Parametric bootstrap likelihood ratio test
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
pbLRT <- function(mod, mod_null, B = 1000) {
  t.observed <- 2 * (mod$logL - mod_null$logL)
  BOOT <- rep(NA, B)
  for (j in 1:B) {
    #print(j)
    BOOT[j] <- doBoot(mod = mod, mod_null = mod_null, test = "LRT")
  }
  perc.rank <- function(x, y) (1 + sum(na.omit(y) >= x)) / (length(na.omit(y)) + 1)
  p.val <- perc.rank(t.observed, BOOT)
  return(p.val)
}
