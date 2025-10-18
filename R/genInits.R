#' Generate initialization for optimization
#'
#' @param W Numeric vector of counts
#' @param M Numeric vector of sequencing depth
#' @param X Matrix of covariates associated with abundance (including intercept)
#' @param X_star Matrix of covariates associated with dispersion (including intercept)
#' @param np Number of covariates associated with abundance (including intercept)
#' @param npstar Number of covariates associated with dispersion (including intercept)
#' @param link ink function for abundance covariates
#' @param phi.link ink function for dispersion covariates
#' @param nstart Integer. Defaults to \code{1}. Number of starts for optimization.
#' @param use Boolean. Defaults to \code{TRUE}. Indicator of whether to use deterministic intialization.
#'
#' @return Matrix of initializations
#'
#' @examples
#' set.seed(1)
#' seq_depth <- rpois(20, lambda = 10000)
#' my_counts <- rbinom(20, size = seq_depth, prob = 0.001) * 10
#' my_covariate <- cbind(rep(c(0,1), each = 10))
#' colnames(my_covariate) <- c("X1")
#'
#' genInits(W = my_counts, M = seq_depth,
#'        X = cbind(1, my_covariate), X_star = cbind(1, my_covariate),
#'        np = 2, npstar = 2,
#'        link = "logit",
#'        phi.link = "logit", nstart = 2, use = TRUE)
#' @export
genInits <- function(W, M,
                     X, X_star,
                     np, npstar,
                     link, phi.link,
                     nstart = 1, use = TRUE) {


  init.glm <- eval(parse(text = paste("quasibinomial(link =", link, ")")))
  tmp <- stats::glm.fit(x = X, y = cbind(W, M - W), family = init.glm)
  b_init <- stats::coef(tmp)
  # Just use 0.5 for phi_init
  # disp_init <- sum((tmp$weights*tmp$residuals^2)[tmp$weights > 0])/tmp$df.r
  # phi_init <- disp_init/(stats::median(M) - 1)

  phi_init <- 0.5
  bstar_int_init <- switch(phi.link, "fishZ" = fishZ(phi_init), "logit" = logit(phi_init))
  bstar_init <- c(bstar_int_init, rep(0, npstar - 1))

  init_start <- rbind(c(b_init, bstar_init))
  if (use) {
    inits <- init_start
  } else {
    inits <- stats::rnorm(length(init_start), init_start, .25)
  }

  if (nstart > 1) {
    for (i in 2:nstart) {
      inits <- rbind(inits, stats::rnorm(length(init_start), init_start, .25))
    }
  }

  return(inits)
}
