#' Maximum Likelihood for the Beta-binomial Distribution
#'
#' @param formula Formula for mean
#' @param phi.formula Formula for overdispersion
#' @param data Data frame
#' @param link Link function for mean, defaults to "logit"
#' @param phi.link Link function for overdispersion, defaults to "fishZ"
#' @param phi.init Optional starting value for overdispersion
#' @param method Method for optimization (see \code{\link{optimr}})
#' @param control Optimization control parameters (see \code{\link{optimr}})
#' @param numerical Indicator to use numerical derivative, useful for testing
#' @param ... Additional arguments for \code{\link{optimr}}
#'
#' @return BBD Model fit
#'
#' @export
bbdml <- function(formula, phi.formula, data,
                  link = "logit",
                  phi.link = "fishZ",
                  phi.init = NULL,
                  method = "L-BFGS-B",
                  control = list(maxit = 1000, reltol = 1e-14, factr = 1e3),
                  numerical = FALSE,
                  ...) {
  if (numerical) {
    control$usenumDeriv <- TRUE
  }

  # Record call
  call <- match.call(expand.dots = FALSE)
  # Record mu link
  link <- match.arg(link)
  # Record phi link
  phi.link <- match.arg(phi.link)
  # Record Method
  method <- match.arg(method)

  mu.f <- formula
  phi.f <- phi.formula

  # If want to model overdispersion, remove intercept for treatment contrast
  # if (phi.f != ~ 1) {
  #   phi.f <- stats::update.formula(phi.f, ~ . - 1)
  # }

  # Get data frame from formula
  dat <- data
  # response
  resp <- stats::model.response(stats::model.frame(formula = mu.f, data = dat))

  # mu model matrix
  X.b <- stats::model.matrix(object = mu.f, data = dat)
  # phi model matrix
  X.bstar <- stats::model.matrix(object = phi.f, data = dat)

  # Number of parameters
  np <- ncol(X.b)
  npstar <- ncol(X.bstar)
  nppar <- np + npstar


  # Get link for mu initialization
  init.glm <- eval(parse(text = paste("binomial(link =", link,")")))
  # Mu initialization
  mu.init.mod <- stats::glm(formula = mu.f, family = init.glm, data = dat)

  if (is.null(phi.init)) {
    z <- 0.1
    # Takes scale, applies inverse to z
    phi.init <- switch(phi.link, "fishZ" = invfishZ(z))
  }

  # Counts
  W <- resp[, 1]
  # Sample Size
  M <- rowSums(resp)

  # Get full initializations
  phi.init <- rep(phi.init, npstar)
  theta.init <- c(stats::coef(mu.init.mod), phi.init)

  # Counts
  W <- resp[, 1]
  # Sample Size
  M <- rowSums(resp)



  if (method == "L-BFGS-B") {
    lower <- c(rep(logit(.0001), np), rep(fishZ(0), npstar))
    upper <- c(rep(logit(.99), np), rep(fishZ(100/(max(M) - 1)), npstar))
    starttime <- proc.time()[1]
    mlout <- try(optimr::optimr(par = theta.init,
                            fn = dbetabin,
                            gr = gr_full,
                            lower = lower,
                            upper = upper,
                            method = method,
                            control = control,
                            W = W,
                            M = M,
                            X = X.b,
                            X_star = X.bstar,
                            np = np,
                            npstar = npstar,
                            logpar = TRUE))
    attempts <- 1
    while (class(mlout) == "try-error" && attempts < 10) {
      theta.init <- theta.init * .95
      mlout <- try(optimr::optimr(par = theta.init,
                                  fn = dbetabin,
                                  gr = gr_full,
                                  lower = lower,
                                  upper = upper,
                                  method = method,
                                  control = control,
                                  W = W,
                                  M = M,
                                  X = X.b,
                                  X_star = X.bstar,
                                  np = np,
                                  npstar = npstar,
                                  logpar = TRUE))
      attempts <- attempts + 1
    }
    if (attempts == 10) {
      stop("Too many initializations!")
    }
    time <- proc.time()[1] - starttime
  }


  ## Results


  theta <- mlout$par
  namb <- colnames(X.b)
  namphi <- paste("phi", colnames(X.bstar), sep = ".")
  names(theta) <- c(namb, namphi)

  b      <- utils::head(theta, np)
  b_star <- utils::tail(theta, npstar)

  # other results
  # if fixpar is not null, df.model is lower than nbpar
  df.model <- length(theta)
  logL <- -mlout$value
  df.residual <- length(M) - df.model
  iterations <- mlout$counts
  code <- mlout$convergence
  msg <- if (!is.null(mlout$message)) mlout$message else character(0)

  structure(
    list(
      call = call, link = link, dat = dat,
      formula = mu.f, phi.formula = phi.f, phi.link = phi.link,
      X.mu = X.b, X.phi = X.bstar,
      resp = resp,
      param = theta, b = b, phi = b_star,
      np = nppar, df.model = df.model, df.residual = df.residual,
      logL = logL,
      iterations = iterations, code = code, msg = msg, time = time),
    class = "bbdml")
}

