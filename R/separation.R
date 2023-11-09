# Wrapper around detectseparation::detect_separation
#
# All arguments are always passed to detect_separation as is, except for x.
# If x has only 1 column, x is replaced with cbind(1, x)
#
# Return value is TRUE or FALSE, extracted from either the "separation" or "outcome" element returned by detect_separation
# The name of this element changed between detectseparation package versions 0.2 and 0.3
separationDetection <- function(y, x, family, control) {
  if (ncol(x) == 1) x <- cbind(1, x)
  separation <- detectseparation::detect_separation(y = y, x = x, family = family, control = control)

  # name of boolean element in list returned by detect_separation changed between package versions
  if (utils::packageVersion("detectseparation") < "0.3") {
    return(separation[["separation"]])
  } else {
    return(separation[["outcome"]])
  }
}

# Simple helper to generate warning about detected separation in a bbdml model
# model_name takes a string, e.g. "abundance model" or "dispersion model"
separationWarning <- function(model_name) {
  warning(paste(
    paste0("Separation detected in ", model_name, "!"),
    "Likely, one of your covariates/experimental conditions is such that",
    "there are all zero counts within a group. The results of this model should",
    "be interpreted with care because there is insufficient data to distinguish between groups. \n",
    sep = "\n"
  ), immediate. = TRUE, call. = FALSE)
}
