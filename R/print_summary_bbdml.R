#' Print summary function
#'
#' @param x Model summary output from \code{\link{summary.bbdml}}
#' @param digits the nnumber of significant digits to use when printing.
#' @param signif.stars logical. If \code{TRUE}, `significance stars' are printed for each coefficient.
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

  cat("\nCoefficients:\n")
  coefs <- x$coefficients
  stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
  cat("\n\nLog-likelihood: ", format(x$logL, digits = max(4L, digits + 1L)), sep = "")
}
