#' Densities of beta binomial distributions, permitting non integer x
#'
#' In some cases we may not have integer W and M's. In these cases,
#' we can still use corncob to estimate parameters, but we need to think of them as
#' no longer coming from the specific beta binomial parametric model,
#' and instead from an estimating equations framework.
#'
#'
#' @author Thomas W Yee
#' @author Xiangjie Xue
#' @author Amy D Willis
#'
#' @param x the value at which defined the density
#' @param size number of trials
#' @param prob the probability of success
#' @param rho the correlation parameter
#' @param log if TRUE, log-densities p are given

dbetabinom_cts <-  function (x, size, prob, rho = 0, log = FALSE) {
  dbetabinom_ab_cts(x = x, size = size, shape1 = prob * (1 - rho)/rho,
                        shape2 = (1 - prob) * (1 - rho)/rho, limit.prob = prob,
                        log = log)
}

dbetabinom_ab_cts <- function (x, size, shape1, shape2, log = FALSE, Inf.shape = exp(20),
                               limit.prob = 0.5) {
  Bigg <- Inf.shape
  Bigg2 <- Inf.shape
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)
  LLL <- max(length(x), length(size), length(shape1), length(shape2),
             length(limit.prob))
  if (length(x) != LLL)
    x <- rep_len(x, LLL)
  if (length(size) != LLL)
    size <- rep_len(size, LLL)
  if (length(shape1) != LLL)
    shape1 <- rep_len(shape1, LLL)
  if (length(shape2) != LLL)
    shape2 <- rep_len(shape2, LLL)
  if (length(limit.prob) != LLL)
    limit.prob <- rep_len(limit.prob, LLL)
  is.infinite.shape1 <- is.infinite(shape1)
  is.infinite.shape2 <- is.infinite(shape2)
  ans <- x
  ans[TRUE] <- log(0)
  ans[is.na(x)] <- NA
  ans[is.nan(x)] <- NaN
  ok0 <- !is.na(shape1) & !is.na(shape2) & !is.na(x) & !is.na(size)

  ### Amy, October 17 2023
  ### MODIFIED FROM
  # okk <- (round(x) == x) & (x >= 0) & (x <= size) & !is.infinite.shape1 &
  #   !is.infinite.shape2 & ok0
  ### TO
  okk <- (x >= 0) & (x <= size) & !is.infinite.shape1 &
    !is.infinite.shape2 & ok0
  ### TO REFLECT ESTIMATING EQUATIONS FOCUS

  if (any(okk)) { ## false, would be true if W was integer
    ans[okk] <- lchoose(size[okk], x[okk]) + lbeta(shape1[okk] +
                                                     x[okk], shape2[okk] + size[okk] - x[okk]) - lbeta(shape1[okk],
                                                                                                       shape2[okk])
    endpt1 <- (x == size) & ((shape1 < 1/Bigg) | (shape2 <
                                                    1/Bigg)) & ok0
    if (any(endpt1)) {
      ans[endpt1] <- lgamma(size[endpt1] + shape1[endpt1]) +
        lgamma(shape1[endpt1] + shape2[endpt1]) - (lgamma(size[endpt1] +
                                                            shape1[endpt1] + shape2[endpt1]) + lgamma(shape1[endpt1]))
    }
    endpt2 <- (x == 0) & ((shape1 < 1/Bigg) | (shape2 < 1/Bigg)) &
      ok0
    if (any(endpt2)) {
      ans[endpt2] <- lgamma(size[endpt2] + shape2[endpt2]) +
        lgamma(shape1[endpt2] + shape2[endpt2]) - (lgamma(size[endpt2] +
                                                            shape1[endpt2] + shape2[endpt2]) + lgamma(shape2[endpt2]))
    }
    endpt3 <- ((Bigg < shape1) | (Bigg < shape2)) & ok0
    if (any(endpt3)) {
      ans[endpt3] <- lchoose(size[endpt3], x[endpt3]) +
        lbeta(shape1[endpt3] + x[endpt3], shape2[endpt3] +
                size[endpt3] - x[endpt3]) - lbeta(shape1[endpt3],
                                                  shape2[endpt3])
    }
  }
  if (!log.arg) {
    ans <- exp(ans)
  }
  ok1 <- !is.infinite.shape1 & is.infinite.shape2
  ok2 <- is.infinite.shape1 & !is.infinite.shape2
  ok3 <- Bigg2 < shape1 & Bigg2 < shape2
  ok4 <- is.infinite.shape1 & is.infinite.shape2
  if (any(ok3, na.rm = TRUE)) { # FALSE
    ok33 <- !is.na(ok3) & ok3
    prob1 <- shape1[ok33]/(shape1[ok33] + shape2[ok33])
    ans[ok33] <- stats::dbinom(x = x[ok33], size = size[ok33], prob = prob1,
                               log = log.arg)
    if (any(ok4)) {
      ans[ok4] <- stats::dbinom(x = x[ok4], size = size[ok4],
                                prob = limit.prob[ok4], log = log.arg)
    }
  }
  if (any(ok1))
    ans[ok1] <- stats::dbinom(x = x[ok1], size = size[ok1], prob = 0,
                              log = log.arg)
  if (any(ok2))
    ans[ok2] <- stats::dbinom(x = x[ok2], size = size[ok2], prob = 1,
                              log = log.arg)
  ans[is.na(shape1) | shape1 < 0] <- NaN
  ans[is.na(shape2) | shape2 < 0] <- NaN
  a.NA <- is.na(x) | is.na(size) | is.na(shape1) | is.na(shape2)
  ans[a.NA] <- NA
  ans
}
