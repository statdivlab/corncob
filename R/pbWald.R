#' Parametric bootstrap Wald test
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null an object of class \code{bbdml}, should be nested within \code{mod}
#' @param B Integer. Defaults to \code{1000}. Number of bootstrap iterations.
#' @param robust Should robust standard errors be used? If not, model-based standard arras are used. Logical, defaults to \code{FALSE}.
#'
#' @return P-value from parametric bootstrap Wald test.
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
#' pbWald(mod1, mod2, B = 50)
#' @export
pbWald <- function(mod, mod_null, B = 1000, robust = FALSE) {

  if (mod$has_noninteger) stop("Can't perform parametric bootstrap with non-integer M or W. ")

  tmp <- getRestrictionTerms(mod = mod, mod_null = mod_null)
  restrictions <- tmp$mu
  restrictions.phi <- tmp$phi
  t.observed <- try(waldchisq_test(mod, restrictions = restrictions, restrictions.phi = restrictions.phi, robust = robust), silent = TRUE)
  if (inherits(t.observed, "try-error")) {
    return(NA)
  }

  BOOT <- rep(NA, B)
  for (j in 1:B) {
    #print(j)
    BOOT[j] <- doBoot(mod = mod, mod_null = mod_null, test = "Wald", robust = robust)
  }
  perc.rank <- function(x, y) (1 + sum(stats::na.omit(y) >= x)) / (length(stats::na.omit(y)) + 1)
  p.val <- perc.rank(t.observed, BOOT)
  return(p.val)
}
