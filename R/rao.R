#' Rao-type chi-squared test (model-based or robust)
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null an object of class \code{bbdml}, should be nested within \code{mod}
#'
#' @return P-value from likelihood ratio test.
#'
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
#' raotest(mod1, mod2)
#' @export
raotest <- function(mod, mod_null) {

  dof.dif <- mod$df.model - mod_null$df.model

  stopifnot(checkNested(mod, mod_null))
  stopifnot(mod$robust == mod_null$robust)
  stopifnot(dof.dif > 0)

  robust <- mod$robust

  #### Don't want to just calculate the likelihood under the null, at the null.
  #### We want need to calculate the likelihood under the alternative, at the null.
  #### Score and hessian only depend on fitted parameters through mu and phi,
  #### so take the alternative model object and just replace mu and phi with fitted values under null.
  ll_full_at_theta0 <- mod
  ll_full_at_theta0$mu.resp <- mod_null$mu.resp
  ll_full_at_theta0$phi.resp <- mod_null$phi.resp

  ll_full_at_theta0$param[names(mod_null$param)] <- mod_null$param
  ll_full_at_theta0$param[setdiff(names(mod$param), names(mod_null$param))] <- 0
  ## assumption is that value under null is zero

  if (robust) {
    inv_fish_info_null <- sand_vcov(ll_full_at_theta0)
  } else {
    if (mod$has_noninteger) {
      warning("Your data has non-integer W or M, and you aren't using robust testing. We will let you do this, but you should consider robust testing instead because your data cannot be assumed to be drawn from a a beta-binomial distribution.")
    }
    inv_fish_info_null <- try(chol2inv(chol(hessian(ll_full_at_theta0, numerical = FALSE))), silent = TRUE)
    if ("try-error" %in% class(inv_fish_info_null)) {
      inv_fish_info_null <- try(solve(hessian(ll_full_at_theta0, numerical = FALSE)), silent = TRUE)
    }
  }

  score_null <- score(ll_full_at_theta0)

  if ("try-error" %in% class(inv_fish_info_null)) {
    warning("Singular Hessian! Cannot calculate p-values in this setting.", immediate. = TRUE)
    pvalue <- NA
  } else {
    chisq.val <- as.numeric(score_null %*% (inv_fish_info_null) %*% score_null)
    pvalue <- stats::pchisq(chisq.val, dof.dif, lower.tail = FALSE)
  }

  return(pvalue)
}


