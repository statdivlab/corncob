#' Full gradient
#'
#' @param ... See \code{\link{gr_beta}} or \code{\link{gr_betast}}
#'
#' @return Gradient of likelihood with respect to all parameters, as vector
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#'
#' @export
gr_full <- function(...) {
  return(c(gr_beta(...), gr_betast(...)))
}
