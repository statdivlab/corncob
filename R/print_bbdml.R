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
  x <- summary(x)
  # Update below with print_summary_bbdml
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nCoefficients associated with abundance:\n")
  coefs.mu <- x$coefficients[1:x$np.mu,]
  rownames(coefs.mu) <- substring(rownames(coefs.mu), 4)
  stats::printCoefmat(coefs.mu, digits = digits, signif.stars = signif.stars,
                      na.print = "NA", ...)

  cat("\n")

  cat("\nCoefficients associated with dispersion:\n")
  coefs.phi <- x$coefficients[(x$np.mu + 1):nrow(x$coefficients),]
  rownames(coefs.phi) <- substring(rownames(coefs.phi), 5)
  stats::printCoefmat(coefs.phi, digits = digits, signif.stars = signif.stars,
                      na.print = "NA", ...)

  cat("\n\nLog-likelihood: ", format(x$logL, digits = max(4L, digits + 1L)), sep = "")
}
