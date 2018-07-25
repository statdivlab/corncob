#' Print ummary function
#'
#' @param object Model summary output from \code{\link{summary.bbdml}}
#' @param ... See details
#'
#' @details ... TODO
#'
#' @return Model summary print
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
print.summary.bbdml <- function(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  } else {
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4L,
                      dimnames = list(cn, colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
  }
  cat("\n\nLog-likelihood: ", format(x$logL, digits = max(4L, digits + 1L)), sep = "")
}
