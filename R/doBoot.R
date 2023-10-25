#' Function to run a bootstrap iteration
#'
#' Internal function. Not intended for users.
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null an object of class \code{bbdml}
#' @param test Character. Hypothesis testing procedure to use. One of \code{"Wald"} or \code{"LRT"} (likelihood ratio test).
#' @param robust Should robust standard errors be used? If not, model-based standard arras are used. Logical, defaults to \code{FALSE}.
#'
#' @return test statistic from one bootstrap iteration
doBoot <- function(mod, mod_null, test, robust = FALSE) {
  # Simulate n samples from the model fit under the null
  newW <- simulate(object = mod_null, nsim = nrow(mod_null$dat))

  newdf <- mod$dat
  newdf$W <- newW
  # Refit models using simulated data
  newout_null <- suppressWarnings(try(bbdml(formula = mod_null$formula,
                                            phi.formula = mod_null$phi.formula,
                                            link = mod_null$link, phi.link = mod_null$phi.link,
                                            data = newdf, inits = mod_null$inits, robust = robust), silent = TRUE))
  newout_alt <- suppressWarnings(try(bbdml(formula = mod$formula,
                                           phi.formula = mod$phi.formula,
                                           link = mod$link, phi.link = mod$phi.link,
                                           data = newdf, inits = mod$inits, robust = robust), silent = TRUE))

  if ("try-error" %in% c(class(newout_null), class(newout_alt))) {
    return(NA)
  }

  # LRT
  if (test == "LRT") {
    if (robust) stop("Amy needs to think about whether robust LRTs are valid!")
    test.stat <- 2 * abs(newout_alt$logL - newout_null$logL)
  } else if (test == "Wald") {
    tmp <- getRestrictionTerms(mod = mod, mod_null = mod_null)
    restrictions <- tmp$mu
    restrictions.phi <- tmp$phi
    test.stat <- try(waldchisq_test(mod = newout_alt, restrictions = restrictions, restrictions.phi = restrictions.phi,
                                    robust = robust), silent = TRUE)
    if (inherits(test.stat, "try-error")) {
      return(NA)
    }
  } else {
    stop(paste("doBoot doesn't know this type of test! 'LRT' and 'Wald' are permitted, but you said", test))
  }
  return(test.stat)
}


