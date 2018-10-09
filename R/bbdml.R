#' Maximum Likelihood for the Beta-binomial Distribution
#'
#' @param formula Formula for mean
#' @param phi.formula Formula for overdispersion
#' @param data Data frame or \code{phyloseq} object
#' @param link Link function for mean, defaults to "logit"
#' @param phi.link Link function for overdispersion, defaults to "logit"
#' @param method Method for optimization (see \code{\link{optimr}})
#' @param control Optimization control parameters (see \code{\link{optimr}})
#' @param numerical Indicator to use numerical derivative, useful for testing
#' @param nstart Number of starts for optimization, defaults to 1
#' @param inits Optional argument to specify initializations, defaults to NULL
#' @param ... Additional arguments for \code{\link{optimr}}
#'
#' @return BBD Model fit
#'
#' @export
bbdml <- function(formula, phi.formula, data,
                  link = "logit",
                  phi.link = "logit",
                  method = "trust",
                  control = list(maxit = 1000, reltol = 1e-14),
                  numerical = FALSE,
                  nstart = 1,
                  inits  = NULL,
                  ...) {
  if (numerical) {
    control$usenumDeriv <- TRUE
  }

  # Convert phyloseq objects
  if ("phyloseq" %in% class(data)) {
    select <- all.vars(formula)[1]
    data <- convert_phylo(data, select = select)
    # Update formula to match convert_phylo specification
    formula <- stats::update(formula, cbind(W, M) ~ .)
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

  # Generate inits
  if (is.null(inits)) {
    inits <- suppressWarnings(genInits(W = W,
                                       M = M,
                                       X = X.b,
                                       X_star = X.bstar,
                                       np = np,
                                       npstar = npstar,
                                       link = link,
                                       phi.link = phi.link))
  } else {
    nstart <- nrow(inits)
    # Or test feasibility of given inits. Same check as in objfun
    for (i in nrow(inits)) {
      # extract matrix of betas (np x 1), first np entries
      b_init      <- utils::head(inits[i,], np)
      # extract matrix of beta stars (npstar x 1), last npstar entries
      b_star_init <- utils::tail(inits[i,], npstar)

      mu.withlink_init <- X.b %*% b_init
      phi.withlink_init <- X.bstar %*% b_star_init
      mu_init <- switch(link, "logit" = invlogit(mu.withlink_init))
      phi_init <- switch(phi.link, "fishZ" = invfishZ(phi.withlink_init), "logit" = invlogit(phi.withlink_init))


      if (any(is.nan(mu_init)) || any(is.nan(phi_init))) {
        warning(paste("Initialization",i,"invalid. Automatically generating new initialization."), immediate. = TRUE)
        inits[i,] <- suppressWarnings(genInits(W = W,
                                               M = M,
                                               X = X.b,
                                               X_star = X.bstar,
                                               np = np,
                                               npstar = npstar,
                                               link = link,
                                               phi.link = phi.link))
      } else {
        val_init <- suppressWarnings(sum(VGAM::dbetabinom(W, M, prob = mu_init, rho = phi_init, log = TRUE)))
        if (is.nan(val_init) || any(phi_init <= sqrt(.Machine$double.eps)) || any(phi_init >= 1 - sqrt(.Machine$double.eps))) {
          warning(paste("Initialization",i,"invalid. Automatically generating new initialization."), immediate. = TRUE)
          inits[i,] <- suppressWarnings(genInits(W = W,
                                                 M = M,
                                                 X = X.b,
                                                 X_star = X.bstar,
                                                 np = np,
                                                 npstar = npstar,
                                                 link = link,
                                                 phi.link = phi.link))
        }
      }

    } ### END FOR: checking inits
  }


  theta.init <- inits[1,]
  if (method == "BFGS") {
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
        curtime <- proc.time()[1] - starttime

        # if the model is improved
        if (mlout$value < bestOut$value) {
          bestOut <- mlout
          time <- curtime
        }
      } ### END IF bfgs
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
      } ### END IF trust
    } ### END FOR - inits
  } ### END IF - nstarts



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



  namb <- paste("mu", colnames(X.b), sep = ".")
  namphi <- paste("phi", colnames(X.bstar), sep = ".")
  names(theta) <- c(namb, namphi)

  b      <- utils::head(theta, np)
  b_star <- utils::tail(theta, npstar)

  mu.withlink <- X.b %*% b
  phi.withlink <- X.bstar %*% b_star
  mu.resp <- switch(link, "logit" = invlogit(mu.withlink))
  phi.resp <- switch(phi.link, "fishZ" = invfishZ(phi.withlink), "logit" = invlogit(phi.withlink))

  theta.resp <- c(mu.resp,phi.resp)

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

