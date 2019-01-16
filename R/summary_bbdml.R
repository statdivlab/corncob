#' Summary function
#'
#' @param object Object of class \code{bbdml}
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
summary.bbdml <- function(object, ...) {
  # For now, Wald test
  coef.table <- waldt(object)
  keep <- match(c("call", "df.model", "df.residual", "logL", "link", "phi.link", "formula", "phi.formula", "np.mu", "np.phi"),
                names(object), 0L)
  ans <- c(object[keep],
           list(coefficients = coef.table))

  class(ans) <- "summary.bbdml"
  return(ans)
}


