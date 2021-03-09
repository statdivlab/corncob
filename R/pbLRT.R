#' Parametric bootstrap likelihood ratio test
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null an object of class \code{bbdml}, should be nested within \code{mod}
#' @param B Integer. Defaults to \code{1000}. Number of bootstrap iterations.
#'
#' @return P-value from parametric bootstrap likelihood ratio test.
#'
#' @examples
#' data(soil_phylum_small)
#' mod1 <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small)
#'
#' mod2 <- bbdml(formula = OTU.1 ~ 1,
#' phi.formula = ~ 1,
#' data = soil_phylum_small)
#' pbLRT(mod1, mod2, B = 100)
#' @export
pbLRT <- function(mod, mod_null, B = 1000) {
  checkNested(mod, mod_null)
  t.observed <- 2 * (mod$logL - mod_null$logL)
  BOOT <- rep(NA, B)
  for (j in 1:B) {
    BOOT[j] <- doBoot(mod = mod, mod_null = mod_null, test = "LRT")
  }
  perc.rank <- function(x, y) (1 + sum(stats::na.omit(y) >= x)) / (length(stats::na.omit(y)) + 1)
  p.val <- perc.rank(t.observed, BOOT)
  return(p.val)
}
