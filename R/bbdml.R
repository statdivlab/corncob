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
                  method = "BFGS",
                  control = list(maxit = 1000, reltol = 1e-14),
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

  # Counts
  W <- resp[, 1]
  # Sample Size
  M <- rowSums(resp)


  # Get link for mu initialization
  #init.glm <- eval(parse(text = paste("binomial(link =", link,")")))
  # Mu initialization
  #mu.init.mod <- stats::glm(formula = mu.f, family = init.glm, data = dat)
  mu.init.mod <- stats::glm.fit(x = X.b, y = resp, family = init.glm)
  # z <- .1
  # mu.init <- switch(link, "logit" = invlogit(z))
  if (is.null(phi.init)) {
    z <- .05
    # Takes scale, applies inverse to z
    phi.init <- switch(phi.link, "fishZ" = invfishZ(z))
  }


  # # Get full initializations
  # if (np > 1) {
  #   #mu.init <- c(mu.init, rep(-10, np - 1))
  #   mu.init <- rep(mu.init, np)
  # }
  if (npstar > 1) {
    #phi.init <- c(phi.init, rep(0, npstar - 1))
    phi.init <- rep(phi.init, npstar)
  }
  #theta.init <- c(mu.init, phi.init)
  theta.init <- c(stats::coef(mu.init.mod),phi.init)


  if (method == "BFGS") {
    #lower <- c(rep(logit(.0001), np), rep(fishZ(0), npstar))
    #upper <- c(rep(logit(.99), np), rep(fishZ(100/(max(M) - 1)), npstar))
    starttime <- proc.time()[1]
    mlout <- optimr::optimr(par = theta.init,
                            fn = dbetabin,
                            gr = gr_full,
                            method = method,
                            control = control,
                            W = W,
                            M = M,
                            X = X.b,
                            X_star = X.bstar,
                            np = np,
                            npstar = npstar,
                            logpar = TRUE)
    theta.orig <- theta.init
    attempts <- 1
    while (mlout$convergence != 0 && attempts < 20) {
      # try going smaller
      theta.init <- theta.init * .95
      mlout <- optimr::optimr(par = theta.init,
                                  fn = dbetabin,
                                  gr = gr_full,
                                  method = method,
                                  control = control,
                                  W = W,
                                  M = M,
                                  X = X.b,
                                  X_star = X.bstar,
                                  np = np,
                                  npstar = npstar,
                                  logpar = TRUE)
      attempts <- attempts + 1
    }
    if (attempts == 20) {
      # reset try going bigger
      attempts <- 1
      theta.init <- theta.orig
      while (mlout$convergence != 0 && attempts < 20) {
        theta.init <- theta.init * 1.05
        mlout <- optimr::optimr(par = theta.init,
                                    fn = dbetabin,
                                    gr = gr_full,
                                    method = method,
                                    control = control,
                                    W = W,
                                    M = M,
                                    X = X.b,
                                    X_star = X.bstar,
                                    np = np,
                                    npstar = npstar,
                                    logpar = TRUE)
        attempts <- attempts + 1
      }
      if (attempts == 20) {
        stop("Too many initializations!")
      }
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

  b.resp <- switch(link, "logit" = invlogit(b))
  b_star.resp <- switch(phi.link, "fishZ" = invfishZ(b_star))
  theta.resp <- c(b.resp,b_star.resp)
  names(theta.resp) <- names(theta)

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
      param = theta, b.mu = b, b.phi = b_star,
      param.response = theta.resp, b.mu.resp = b.resp, b.phi.resp = b_star.resp,
      np.total = nppar, np.mu = np, np.phi = npstar,
      df.model = df.model, df.residual = df.residual,
      logL = logL,
      iterations = iterations, code = code, msg = msg, time = time),
    class = "bbdml")
}

