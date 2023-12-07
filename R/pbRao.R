#' Parametric bootstrap Rao test
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null an object of class \code{bbdml}, should be nested within \code{mod}
#' @param B Integer. Defaults to \code{1000}. Number of bootstrap iterations.
#'
#' @return P-value from parametric bootstrap Rao test.
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
#' pbRao(mod1, mod2, B = 10)
#' @export
pbRao <- function(mod, mod_null, B = 1000) {

  if (mod$has_noninteger) stop("Can't perform parametric bootstrap with non-integer M or W. ")

  checkNested(mod, mod_null)

  t.observed <- stats::qchisq(raotest(mod, mod_null),
                              mod$df.model - mod_null$df.model,
                              lower.tail = FALSE)
  BOOT <- rep(NA, B)
  for (j in 1:B) {
    BOOT[j] <- doBoot(mod = mod, mod_null = mod_null, test = "Rao")
  }
  perc.rank <- function(x, y) (1 + sum(stats::na.omit(y) >= x)) / (length(stats::na.omit(y)) + 1)
  p.val <- perc.rank(t.observed, BOOT)
  return(p.val)
}
