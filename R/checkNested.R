#' Check for nested models
#'
#' @param mod model fit from bbdml
#' @param mod_null model fit from bbdml
#'
#' @return Boolean. \code{TRUE} if \code{mod_null} is nested within \code{mod}.
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#' @export
checkNested <- function(mod, mod_null) {

  if (!all(names(mod_null$b.mu) %in% names(mod$b.mu))) {
    stop("Models for abundance are not nested!")
  }

  if (!all(names(mod_null$b.phi) %in% names(mod$b.phi))) {
    stop("Models for dispersion are not nested!")
  }

  if (!all(mod$dat == mod_null$dat)) {
    stop("Models are not fit to the same data!")
  }

  return(TRUE)
}
