% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbWald.R
\name{pbWald}
\alias{pbWald}
\title{Parametric bootstrap Wald test}
\usage{
pbWald(mod, mod_null, B = 1000, robust = FALSE)
}
\arguments{
\item{mod}{an object of class \code{bbdml}}

\item{mod_null}{an object of class \code{bbdml}, should be nested within \code{mod}}

\item{B}{Integer. Defaults to \code{1000}. Number of bootstrap iterations.}

\item{robust}{Should robust standard errors be used? If not, model-based standard arras are used. Logical, defaults to \code{FALSE}.}
}
\value{
P-value from parametric bootstrap Wald test.
}
\description{
Parametric bootstrap Wald test
}
\examples{
data(soil_phylum_small_otu1)
mod1 <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)

mod2 <- bbdml(formula = cbind(W, M - W) ~ 1,
phi.formula = ~ 1,
data = soil_phylum_small_otu1)
pbWald(mod1, mod2, B = 50)
}
