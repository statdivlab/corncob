#' Summary function
#'
#' @param object Model fit output from \code{\link{bbdml}}
#' @param ... See details
#'
#' @details ... TODO
#'
#' @return Model summary
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
summary.bbdml <- function(object, ...) {
  # For now, just Wald test
  coef.table <- waldtest(object)
  keep <- match(c("call", "df.model", "df.residual", "logL", "link", "phi.link", "formula", "phi.formula"),
                names(object), 0L)
  ans <- c(object[keep],
           list(coefficients = coef.table,
                aliased = aliased))

  class(ans) <- "summary.bbdml"
  return(ans)
}


