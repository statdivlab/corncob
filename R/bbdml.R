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
#' @param nstart Number of starts for optimization, defaults to 10
#' @param inits Optional argument to specify initializations, defaults to NULL
#' @param ... Additional arguments for \code{\link{optimr}}
#'
#' @return BBD Model fit
#'
#' @export
bbdml <- function(formula, phi.formula, data,
                  link = "logit",
                  phi.link = "fishZ",
                  phi.init = NULL,
                  method = "trust",
                  control = list(maxit = 1000, reltol = 1e-14),
                  numerical = FALSE,
                  nstart = 10,
                  inits  = NULL,
                  ...) {
  if (numerical) {
    control$usenumDeriv <- TRUE
  }

  # Record call
  call <- match.call(expand.dots = FALSE)
  # Record mu link
  link <- match.arg(link, choices = "logit")
  # Record phi link
  phi.link <- match.arg(phi.link, choices = c("fishZ", "logit"))
  # Record Method
  method <- match.arg(method, choices = c("BFGS", "trust"))

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
#### COMMENTED BELOW FOR NEW INIT
#
#   # Get link for mu initialization
#   init.glm <- eval(parse(text = paste("binomial(link =", link,")")))
#   # Mu initialization
#   #mu.init.mod <- stats::glm(formula = mu.f, family = init.glm, data = dat)
#   # Add fake counts no just for stable initializations
#   fakeresp <- resp
#   zhold <- which(fakeresp[,1] == 0)
#   fakeresp[zhold, 1] <- 1
#   fakeresp[zhold, 2] <- fakeresp[zhold, 2] - 1
#   init.mod <- stats::glm.fit(x = X.b, y = fakeresp, family = init.glm)
#   mu.init <- stats::coef(init.mod)
#   # fix for numerical stability
#   if (any(invlogit(mu.init) == 0)) {
#     mu.init[invlogit(mu.init) == 0] <- logit(1/mean(M))
#   }
#   if (any(invlogit(mu.init) == 1)) {
#     mu.init[invlogit(mu.init) == 1] <- logit(1 - 1/mean(M))
#   }
#   # z <- .1
#   # mu.init <- switch(link, "logit" = invlogit(z))
#   # won't mess with links because should be close to 0
#   phi.init <- 1/(mean(M) - 1)
#
#   # # Get full initializations
#   # if (np > 1) {
#   #   #mu.init <- c(mu.init, rep(-10, np - 1))
#   #   mu.init <- rep(mu.init, np)
#   # }
#   if (npstar > 1) {
#     #phi.init <- c(phi.init, rep(0, npstar - 1))
#     phi.init <- c(phi.init, rep(0, npstar - 1))
#   }
#   #theta.init <- c(mu.init, phi.init)
#   theta.init <- c(mu.init, phi.init)
  if (is.null(inits)) {
    inits <- suppressWarnings(genInits(nstart = nstart,
                                       W = W,
                                       M = M,
                                       X = X.b,
                                       X_star = X.bstar,
                                       np = np,
                                       npstar = npstar,
                                       link = link,
                                       phi.link = phi.link,
                                       logpar = TRUE))
  }


  theta.init <- inits[1,]
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
                            link = link,
                            phi.link = phi.link,
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
                                  link = link,
                                  phi.link = phi.link,
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
                                    link = link,
                                    phi.link = phi.link,
                                    logpar = TRUE)
        attempts <- attempts + 1
      }
      if (attempts == 20) {
        stop("Too many initializations!")
      }
    }
    curtime <- proc.time()[1] - starttime
  }
  if (method == "trust") {
    if (!exists("rinit")) {
      rinit <- 1
    }
    if (!exists("rmax")) {
      rmax <- 100
    }
    starttime <- proc.time()[1]
    mlout <- trust::trust(objfun, parinit = theta.init,
                          W = W,
                          M = M,
                          X = X.b,
                          X_star = X.bstar,
                          np = np,
                          npstar = npstar,
                          link = link,
                          phi.link = phi.link,
                          rinit = rinit,
                          rmax = rmax)
    curtime <- proc.time()[1] - starttime
  }
  # Save the best model
  bestOut <- mlout
  time <- curtime

  if (nstart >= 2) {
    for (i in 2:nstart) {
      ### BEGIN FOR
      theta.init <- inits[i,]
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
                                link = link,
                                phi.link = phi.link,
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
                                  link = link,
                                  phi.link = phi.link,
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
                                    link = link,
                                    phi.link = phi.link,
                                    logpar = TRUE)
            attempts <- attempts + 1
          }
          if (attempts == 20) {
            stop("Too many initializations!")
          }
        }
        curtime <- proc.time()[1] - starttime

        # if the model is improved
        if (mlout$value < bestOut$value) {
          bestOut <- mlout
          time <- curtime
        }
      }
      if (method == "trust") {
        starttime <- proc.time()[1]
        mlout <- trust::trust(objfun, parinit = theta.init,
                              W = W,
                              M = M,
                              X = X.b,
                              X_star = X.bstar,
                              np = np,
                              npstar = npstar,
                              link = link,
                              phi.link = phi.link,
                              rinit = rinit,
                              rmax = rmax)
        curtime <- proc.time()[1] - starttime

        # if the model is improved
        if (mlout$value < bestOut$value) {
          bestOut <- mlout
          time <- curtime
        }
      }


      ### END FOR - inits
    }
  }



  # change back for name
  mlout <- bestOut

  if (method == "BFGS") {
    theta <- mlout$par
    logL <- -mlout$value
    iterations <- mlout$counts
    code <- mlout$convergence
  }
  if (method == "trust") {
    theta <- mlout$argument
    logL <- -mlout$value
    iterations <- mlout$iterations
    code <- mlout$converged
  }
  ## Results



  namb <- colnames(X.b)
  namphi <- paste("phi", colnames(X.bstar), sep = ".")
  names(theta) <- c(namb, namphi)

  b      <- utils::head(theta, np)
  b_star <- utils::tail(theta, npstar)

  mu.withlink <- X.b %*% b
  phi.withlink <- X.bstar %*% b_star
  mu.resp <- switch(link, "logit" = invlogit(mu.withlink))
  phi.resp <- switch(phi.link, "fishZ" = invfishZ(phi.withlink), "logit" = invlogit(phi.withlink))

  theta.resp <- c(mu.resp,phi.resp)
  #names(theta.resp) <- names(theta)

  # other results
  # if fixpar is not null, df.model is lower than nbpar
  df.model <- length(theta)

  df.residual <- length(M) - df.model

  msg <- if (!is.null(mlout$message)) mlout$message else character(0)

  structure(
    list(
      call = call, link = link, dat = dat,
      formula = mu.f, phi.formula = phi.f, phi.link = phi.link,
      X.mu = X.b, X.phi = X.bstar,
      resp = resp, M = M, W = W,
      param = theta, b.mu = b, b.phi = b_star,
      param.response = theta.resp, mu.resp = mu.resp, phi.resp = phi.resp,
      np.total = nppar, np.mu = np, np.phi = npstar,
      df.model = df.model, df.residual = df.residual,
      logL = logL, inits = inits,
      iterations = iterations, code = code, msg = msg, time = time),
    class = "bbdml")
}

