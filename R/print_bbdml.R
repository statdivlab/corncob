#' Print function
#'
#' @param x Object of class \code{bbdml}
#' @param digits the nnumber of significant digits to use when printing.
#' @param signif.stars logical. If \code{TRUE}, `significance stars' are printed for each coefficient.
#' @param ... No optional arguments are accepted at this time.
#'
# #' @details ...
#'
#' @return Model summary
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
print.bbdml <- function(x, digits = max(3L, getOption("digits") - 3L),
                        signif.stars = getOption("show.signif.stars"), ...) {
  # Update with print_summary_bbdml
  x <- summary(x)
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nCoefficients:\n")
  coefs <- x$coefficients
  stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                      na.print = "NA", ...)
  cat("\n\nLog-likelihood: ", format(x$logL, digits = max(4L, digits + 1L)), sep = "")
}
