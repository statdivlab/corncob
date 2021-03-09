#' Summary function
#'
#' @param object Object of class \code{bbdml}
#' @param ... No optional arguments are accepted at this time.
#'
#'
#' @return Object of class \code{summary.bbdml}. Displays printed model summary.
#'
#' @examples
#' data(soil_phylum_small)
#' mod <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil_phylum_small)
#' summary(mod)
#' @export
summary.bbdml <- function(object, ...) {
  # For now, Wald test
  coef.table <- waldt(object)
  # keep <- match(c("call", "df.model", "df.residual", "logL", "link", "phi.link", "formula", "phi.formula", "np.mu", "np.phi", "sep_da", "sep_dv"),
  #               names(object), 0L)
  keep <- match(c("call", "df.model", "df.residual", "logL", "link", "phi.link", "formula", "phi.formula", "np.mu", "np.phi", "sep_da", "sep_dv", "mu.resp", "phi.resp"),
                names(object), 0L)
  ans <- c(object[keep],
           list(coefficients = coef.table))

  class(ans) <- "summary.bbdml"
  return(ans)
}


