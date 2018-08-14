#' Generate initialization for optimization
#'
#' @param W absolute abundance
#' @param M sample size
#' @param X mean covariates
#' @param X_star overdispersion covariates
#' @param np number of mean parameters
#' @param npstar number of overdisperion parameters
#' @param link Link function for mean
#' @param phi.link Link function for overdispersion
#'
#' @return Matrix of initializations
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
genInits <- function(W, M,
                     X, X_star,
                     np, npstar,
                     link, phi.link) {


  init.glm <- eval(parse(text = paste("quasibinomial(link =", link,")")))
  tmp <- stats::glm.fit(x = X, y = cbind(W, M - W), family = init.glm)
  b_init <- stats::coef(tmp)
  disp_init <- sum((tmp$weights*tmp$residuals^2)[tmp$weights > 0])/tmp$df.r
  phi_init <- disp_init/(stats::median(M) - 1)
  bstar_int_init <- switch(phi.link, "fishZ" = fishZ(phi_init), "logit" = logit(phi_init))
  bstar_init <- c(bstar_int_init, rep(0, npstar - 1))
  inits <- rbind(c(b_init, bstar_init))

  # if (is.null(lower)) {
  #   lower <- rep(-20, np + npstar)
  # }
  # if (is.null(upper)) {
  #   upper <- rep(20, np + npstar)
  # }
  # if (is.null(max.call)) {
  #   max.call <- 1000
  # }
  # if (is.null(temperature)) {
  #   temperature <- 50000
  # }
#  inits <- matrix(NA, nrow = nstart, ncol = np + npstar)
#   for (i in 1:nstart) {
#     inits[i,] <- try(GenSA::GenSA(fn = dbetabin, lower = lower, upper = upper, W = W, M = M, X = X, X_star = X_star, np = np,
#                          npstar = npstar, link = link, phi.link = phi.link,
#                          control = list(max.call = max.call, temperature = temperature, smooth = FALSE)), silent = TRUE)$par
#
#     # Add check, same as in objfun and hessian
#     # extract matrix of betas (np x 1), first np entries
#     b_init      <- utils::head(inits[i,], np)
#     # extract matrix of beta stars (npstar x 1), last npstar entries
#     b_star_init <- utils::tail(inits[i,], npstar)
#
#     mu.withlink_init <- X %*% b_init
#     phi.withlink_init <- X_star %*% b_star_init
#     mu_init <- switch(link, "logit" = invlogit(mu.withlink_init))
#     phi_init <- switch(phi.link, "fishZ" = invfishZ(phi.withlink_init), "logit" = invlogit(phi.withlink_init))
#
#     val_init <- suppressWarnings(sum(VGAM::dbetabinom(W, M, prob = mu_init, rho = phi_init, log = TRUE)))
#     if (is.nan(val_init) || any(phi_init <= sqrt(.Machine$double.eps)) || any(phi_init >= 1 - sqrt(.Machine$double.eps))) {
#       stop("Cannot generate initializations! \n\n You are likely overparameterizing phi.formula \n without enough signal in the data. \n\n Consider setting phi.link = \"logit\",\n removing covariates in phi.formula,\n or setting your own initializations. \n\n If none of the above works, you probably \n have underdispersion in your data. \n This cannot be modeled with a beta-binomial. \n Consider using a binomial GLM, or use \n a quasibinomial GLM to explicitly \n model underdispersion.")
#     }
# }

  return(inits)
}
