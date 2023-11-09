#' Print summary function
#'
#' @param x Object of class \code{bbdml}
#' @param digits the number of significant digits to use when printing.
#' @param signif.stars logical. If \code{TRUE}, `significance stars' are printed for each coefficient.
#' @param ... No optional arguments are accepted at this time.
#'
#'
#' @return \code{NULL}. Displays printed model summary.
#'
#' @examples
#' data(soil_phylum_small)
#' mod <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small)
#' print(summary(mod))
#' @export
print.summary.bbdml <- function(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nCoefficients associated with abundance:\n")
  coefs.mu <- x$coefficients[1:x$np.mu,, drop = FALSE]
  rownames(coefs.mu) <- substring(rownames(coefs.mu), 4)
  stats::printCoefmat(coefs.mu, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)

  cat("\n")

  cat("\nCoefficients associated with dispersion:\n")
  coefs.phi <- x$coefficients[(x$np.mu + 1):nrow(x$coefficients),, drop = FALSE]
  rownames(coefs.phi) <- substring(rownames(coefs.phi), 5)
  stats::printCoefmat(coefs.phi, digits = digits, signif.stars = signif.stars,
                      na.print = "NA", ...)

  cat("\n\nLog-likelihood: ", format(x$logL, digits = max(4L, digits + 1L)),
      "\n", sep = "")

  if (x$sep_da || x$sep_dv) {
    warning("This model is based on a discriminant taxa.
You may see NAs in the model summary because Wald testing is invalid.
Likelihood ratio and Rao testing can be used, but valid standard errors cannot be calculated.",
            immediate. = TRUE)
  }
}
