#' Logit transformation
#'
#' @param x data
#'
#' @return logit of x
#'
#' @examples
#' x <- .5
#' logit(x)
#'
#' @export
logit <- function(x) {
  return(stats::qlogis(x))
}


#' Fisher's z transformation
#'
#' @param x data
#'
#' @return Fisher's z transformation of x
#'
#' @examples
#' x <- .5
#' fishZ(x)
#'
#' @export
fishZ <- function(x) {
  #return(log((1 + x)/(1 - x)))
  return(2 * atanh(x))
}

#' Inverse logit transformation
#'
#' @param x data
#'
#' @return Inverse logit transformation of x
#'
#' @examples
#' x <- .5
#' invlogit(x)
#'
#' @export
invlogit <- function(x) {
  return(stats::plogis(x))
}

#' Inverse Fisher's z transformation
#'
#' @param x data
#'
#' @return Inverse Fisher's z transformation of x
#'
#' @examples
#' x <- .5
#' invfishZ(x)
#'
#' @export
invfishZ <- function(x) {
  #return((exp( x) - 1)/(1 + exp(x)))
  return(tanh(0.5 * x))
}

#' Hyperbolic cotangent transformation
#'
#' @param x data
#'
#' @return Hyperbolic cotangent transformation of x
#'
#' @examples
#' x <- .5
#' coth(x)
#'
#' @export
coth <- function(x) {
  return(c(1/tanh(x)))
}
