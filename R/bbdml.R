#' Maximum Likelihood for the Beta-binomial Distribution
#'
#' @param formula an object of class \code{formula}: a symbolic description of the model to be fitted to the abundance
#' @param phi.formula an object of class \code{formula} without the response: a symbolic description of the model to be fitted to the dispersion
#' @param data a data frame or \code{phyloseq} object containing the variables in the models
#' @param link link function for abundance covariates, defaults to \code{"logit"}
#' @param phi.link link function for dispersion covariates, defaults to \code{"logit"}
#' @param method optimization method, defaults to \code{"trust"}, or see \code{\link[optimx]{optimr}} for other options
#' @param control optimization control parameters (see \code{\link[optimx]{optimr}})
#' @param numerical Boolean. Defaults to \code{FALSE}. Indicator of whether to use the numeric Hessian (not recommended).
#' @param nstart Integer. Defaults to \code{1}. Number of starts for optimization.
#' @param inits Optional initializations as rows of a matrix. Defaults to \code{NULL}.
#' @param allow_noninteger Boolean. Defaults to \code{FALSE}. Should noninteger W's and M's be allowed? This behavior was not permitted prior to v4.1, needs to be explicitly allowed.
#' @param robust Should robust standard errors be returned? If not, model-based standard arras are used. Logical, defaults to \code{FALSE}.
#' @param ... Optional additional arguments for \code{\link[optimx]{optimr}} or \code{\link{trust}}
#'
#' @return An object of class \code{bbdml}.
#'
#' @examples
#' # data frame example
#' data(soil_phylum_small_otu1)
#' bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small_otu1)
#'
#' # phyloseq example (only run this if you have phyloseq installed)
#' \dontrun{
#' data(soil_phylum_small_sample)
#' data(soil_phylum_small_otu)
#' data_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylum_small_sample),
#' phyloseq::otu_table(soil_phylum_small_otu, taxa_are_rows = TRUE))
#' bbdml(formula = Proteobacteria ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = data_phylo)
#' }
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
                  allow_noninteger = FALSE,
                  robust = FALSE,
                  ...) {
  if (numerical) {
    control$usenumDeriv <- TRUE
  }

  argList <- list(...)

  # Convert phyloseq objects
  if ("phyloseq" %in% class(data)) {
    selection <- all.vars(formula)[1]
    data <- convert_phylo(data, select = selection)
    # Update formula to match convert_phylo specification
    formula <- stats::update(formula, cbind(W, M - W) ~ .)
  }

  # Record call
  call <- match.call(expand.dots = FALSE)
  # Record mu link
  link <- try(match.arg(link, choices = "logit"))
  if ("try-error" %in% class(link)) {
    stop('link currently only allows "logit"! \n
         To request another link function, post an issue on GitHub.')
  }

  # Record phi link
  phi.link <- try(match.arg(phi.link, choices = c("fishZ", "logit")))
  if ("try-error" %in% class(phi.link)) {
    stop('phi.link currently only allows "logit" or "fishZ"! \n
         To request another link function, post an issue on GitHub.')
  }

  # Record Method
  method <- try(match.arg(method, choices = c("BFGS", "trust")))
  if ("try-error" %in% class(method)) {
    stop('If method is specified, it must be either "BFGS" or "trust"!')
  }
  # Check that optimx is installed if method is "BFGS"
  if (method == "BFGS") {
    optimx_install <- try(find.package("optimx"), silent = TRUE)
    if (inherits(optimx_install, "try-error")) {
      {stop("If you would like to use the 'BFGS' method, please install the `optimx` package from CRAN.")}
    }
  }

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

  # mu terms
  terms.mu <- stats::terms(mu.f)
  # phi terms
  terms.phi <- stats::terms(phi.f)

  # Number of parameters
  np <- ncol(X.b)
  npstar <- ncol(X.bstar)
  nppar <- np + npstar

  df.model <- nppar
  df.residual <- nrow(resp) - df.model

  if (df.residual <= 0) {
    stop("Model overspecified! \n
Trying to fit more parameters than sample size. Model cannot be estimated.")
  }


  # Counts
  W <- resp[, 1]
  # Sample Size
  M <- rowSums(resp)

  has_noninteger <- any(round(W) != W) | any(round(M) != M)
  if (!allow_noninteger & has_noninteger) {
    stop("You don't have counts in your abundances. \n
         Double-check and correct your data, and/or force the model to fit with `allow_noninteger`")
  }

  # Check for separation
  sep_da <- sep_dv <- FALSE
  if (length(attr(terms.mu, "term.labels") != 0)) {

    withCallingHandlers(
      expr = {
        sep_da <- separationDetection(
          y = cbind(W, M - W), x = X.b, family = stats::binomial("logit"),
          control = list(implementation = "lpSolveAPI",
                         linear_program = "primal",
                         purpose = "test")
        )
      },
      warning = function(w) {
        if (allow_noninteger & grepl("non-integer counts in a binomial glm", w$message)) {
          invokeRestart("muffleWarning")
        } else {
          warning(w$message, call.=F)
          invokeRestart("muffleWarning")
        }
      }
    )

    if (sep_da) separationWarning(model_name = "abundance model")
  }

  if (length(attr(terms.phi, "term.labels") != 0)) {

    withCallingHandlers(
        expr = {
          sep_dv <- separationDetection(
            y = cbind(W, M - W), x = X.bstar, family = stats::binomial("logit"),
            control = list(implementation = "lpSolveAPI",
                           linear_program = "primal",
                           purpose = "test")
          )
        },
        warning = function(w) {
          if (allow_noninteger & grepl("non-integer counts in a binomial glm", w$message)) {
            invokeRestart("muffleWarning")
          } else {
            warning(w$message, call.=F)
            invokeRestart("muffleWarning")
          }
        }
      )

    if (sep_dv) separationWarning(model_name = "dispersion model")
  }

  # Generate inits
  if (is.null(inits)) {
    inits <- suppressWarnings(genInits(W = W,
                                       M = M,
                                       X = X.b,
                                       X_star = X.bstar,
                                       np = np,
                                       npstar = npstar,
                                       link = link,
                                       phi.link = phi.link,
                                       nstart = nstart,
                                       use = TRUE))
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
                                               phi.link = phi.link,
                                               nstart = 1,
                                               use = FALSE))
      } else {
        val_init <- suppressWarnings(sum(dbetabinom_cts(W, M, prob = mu_init, rho = phi_init, log = TRUE)))
        if (is.nan(val_init) || any(phi_init <= sqrt(.Machine$double.eps)) || any(phi_init >= 1 - sqrt(.Machine$double.eps))) {
          warning(paste("Initialization",i,"invalid. Automatically generating new initialization."), immediate. = TRUE)
          inits[i,] <- suppressWarnings(genInits(W = W,
                                                 M = M,
                                                 X = X.b,
                                                 X_star = X.bstar,
                                                 np = np,
                                                 npstar = npstar,
                                                 link = link,
                                                 phi.link = phi.link,
                                                 nstart = 1,
                                                 use = FALSE))
        }
      }

    } ### END FOR: checking inits
  }



  # theta.init <- inits[1,]
  # # replace any NA inits with 0
  # theta.init[which(is.na(theta.init))] <- 0
  # if (method == "BFGS") {
  #   #starttime <- proc.time()[1]
  #   mlout <- optimr::optimr(par = theta.init,
  #                           fn = dbetabin_neg,
  #                           gr = gr_full,
  #                           method = method,
  #                           control = control,
  #                           W = W,
  #                           M = M,
  #                           X = X.b,
  #                           X_star = X.bstar,
  #                           np = np,
  #                           npstar = npstar,
  #                           link = link,
  #                           phi.link = phi.link,
  #                           logpar = TRUE)
  #   theta.orig <- theta.init
  #   #curtime <- proc.time()[1] - starttime
  # }
  # if (method == "trust") {
  #   if (!exists("rinit")) {
  #     rinit <- 1
  #   }
  #   if (!exists("rmax")) {
  #     rmax <- 100
  #   }
  #   #starttime <- proc.time()[1]
  #   mlout <- trust::trust(objfun, parinit = theta.init,
  #                         W = W,
  #                         M = M,
  #                         X = X.b,
  #                         X_star = X.bstar,
  #                         np = np,
  #                         npstar = npstar,
  #                         link = link,
  #                         phi.link = phi.link,
  #                         rinit = rinit,
  #                         rmax = rmax)
  #   #curtime <- proc.time()[1] - starttime
  # }
  # # Save the best model
  # bestOut <- mlout
  #time <- curtime

  bestOut <- NULL

  # if (nstart >= 2) {
  for (i in 1:nstart) {
    ### BEGIN FOR
    theta.init <- inits[i,]
    # replace any NA inits with 0
    theta.init[which(is.na(theta.init))] <- 0
    if (method == "BFGS") {
      #starttime <- proc.time()[1]
      # optimx::optimr prints out all of the control parameters, suppress these print statements with `invisible(capture.output())`
      invisible(utils::capture.output(mlout <- suppressWarnings(try(optimx::optimr(par = theta.init,
                                                                                   fn = dbetabin_neg,
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
                                                                                   phi.link = phi.link, logpar = TRUE),
                                                                    silent = TRUE)))); if (inherits(mlout, "try-error")) next
      #theta.orig <- theta.init
      #curtime <- proc.time()[1] - starttime

      if (is.null(bestOut)) bestOut <- mlout
      # if the model is improved
      if (mlout$value < bestOut$value) {
        bestOut <- mlout
        #time <- curtime
      }
    } ### END IF bfgs
    if (method == "trust") {
      #starttime <- proc.time()[1]
      if (is.null(argList$rinit)) {
        rinit <- 1
      }
      if (is.null(argList$rmax)) {
        rmax <- 100
      }
      mlout <- try(trust::trust(objfun, parinit = theta.init,
                                W = W,
                                M = M,
                                X = X.b,
                                X_star = X.bstar,
                                np = np,
                                npstar = npstar,
                                link = link,
                                phi.link = phi.link,
                                rinit = rinit, rmax = rmax),
                   silent = TRUE); if (inherits(mlout, "try-error")) next
      #curtime <- proc.time()[1] - starttime

      if (is.null(bestOut)) bestOut <- mlout
      # if the model is improved
      if (mlout$value < bestOut$value) {
        bestOut <- mlout
        #time <- curtime
      }
    } ### END IF trust
  } ### END FOR - inits
  # } ### END IF - nstarts

  if (is.null(bestOut)) stop("Model could not be optimized! Try changing initializations or simplifying your model.")

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


  msg <- if (!is.null(mlout$message)) mlout$message else character(0)

  structure(
    list(
      call = call, link = link, dat = dat,
      formula = mu.f, phi.formula = phi.f, phi.link = phi.link,
      X.mu = X.b, X.phi = X.bstar,
      resp = resp, M = M, W = W,
      param = theta, b.mu = b, b.phi = b_star,
      param.response = theta.resp, terms.mu = terms.mu, terms.phi = terms.phi,
      mu.resp = mu.resp, phi.resp = phi.resp,
      np.total = nppar, np.mu = np, np.phi = npstar,
      df.model = df.model, df.residual = df.residual,
      logL = logL, inits = inits, sep_da = sep_da, sep_dv = sep_dv,
      iterations = iterations, code = code, msg = msg,
      robust = robust,
      allow_noninteger = allow_noninteger,
      has_noninteger = has_noninteger),
    class = "bbdml")
}


