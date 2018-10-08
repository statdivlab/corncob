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

  # Default to simpler model if no observed counts
  if (sum(newW) == 0) {
    return(0)
  }
  newdf <- mod$dat
  newdf$W <- newW
  # Refit models using simulated data
  newout_null <- bbdml(formula = mod_null$formula,
                       phi.formula = mod_null$phi.formula,
                       link = mod_null$link, phi.link = mod_null$phi.link,
                       data = newdf)
  newout_alt <- bbdml(formula = mod$formula,
                      phi.formula = mod$phi.formula,
                      link = mod$link, phi.link = mod$phi.link,
                      data = newdf)

  # LRT
  if (test == "LRT") {
    test.stat <- 2 * abs(newout_alt$logL - newout_null$logL)
  } else if (test == "Wald") {
    altterms.mu <- attr(stats::terms(stats::model.frame(newout_alt$formula, data = newout_alt$dat)), "term.labels")
    altterms.phi <- attr(stats::terms(stats::model.frame(newout_alt$phi.formula, data = newout_alt$dat)), "term.labels")
    nullterms.mu <- attr(stats::terms(stats::model.frame(newout_null$formula, data = newout_null$dat)), "term.labels")
    nullterms.phi <- attr(stats::terms(stats::model.frame(newout_null$phi.formula, data = newout_null$dat)), "term.labels")

    res_mu <- setdiff(altterms.mu, nullterms.mu)
    res_phi <- setdiff(altterms.phi, nullterms.phi)
    if (length(res_mu) == 0) {
      testonly <- "phi"
    } else if (length(res_phi) == 0) {
      testonly <- "mu"
    } else {
      testonly <- NULL
    }
    restrictions <- union(res_mu, res_phi)
    test.stat <- waldchisq_test(mod = newout_alt, restrictions = restrictions, testonly = testonly)
  }
  return(test.stat)
}


