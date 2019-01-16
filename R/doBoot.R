#' Function to run a bootstrap iteration
#'
#' @param mod Unrestricted model fit from \code{bbdml}.
#' @param mod_null Restricted model fit from \code{bbdml}.
#' @param test Character string. Type of test, either "LRT" or "Wald".
#'
#' @return test statistic from one bootstrap iteration
doBoot <- function(mod, mod_null, test) {
  # Simulate n samples from the model fit under the null
  newW <- simulate(object = mod_null, nsim = nrow(mod_null$dat))

  newdf <- mod$dat
  newdf$W <- newW
  # Refit models using simulated data
  newout_null <- try(bbdml(formula = mod_null$formula,
                       phi.formula = mod_null$phi.formula,
                       link = mod_null$link, phi.link = mod_null$phi.link,
                       data = newdf, inits = mod_null$inits), silent = TRUE)
  newout_alt <- try(bbdml(formula = mod$formula,
                      phi.formula = mod$phi.formula,
                      link = mod$link, phi.link = mod$phi.link,
                      data = newdf, inits = mod$inits), silent = TRUE)

  if ("try-error" %in% c(class(newout_null), class(newout_alt))) {
    return(NA)
  }

  # LRT
  if (test == "LRT") {
    test.stat <- 2 * abs(newout_alt$logL - newout_null$logL)
  } else if (test == "Wald") {
    tmp <- getRestrictionTerms(mod = mod, mod_null = mod_null)
    restrictions <- tmp$mu
    restrictions.phi <- tmp$phi
    test.stat <- try(waldchisq_test(mod = newout_alt, restrictions = restrictions, restrictions.phi = restrictions.phi), silent = TRUE)
    if (class(test.stat) == "try-error") {
      return(NA)
    }
  }
  return(test.stat)
}


