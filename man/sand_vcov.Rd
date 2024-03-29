% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sandSE.R
\name{sand_vcov}
\alias{sand_vcov}
\title{Compute sandwich estimate of variance-covariance matrix}
\usage{
sand_vcov(mod, numerical = FALSE)
}
\arguments{
\item{mod}{an object of class \code{bbdml}}

\item{numerical}{Boolean. Defaults to \code{FALSE}. Indicator of whether to use the numeric Hessian and score (not recommended).}
}
\value{
Sandwich variance-covariance matrix. \eqn{\hat{A}^{-1} \hat{B} \hat{A}^{-1}}.
}
\description{
Compute sandwich estimate of variance-covariance matrix
}
\examples{
data(soil_phylum_small_otu1)
mod <- bbdml(formula = cbind(W, M - W) ~ DayAmdmt,
phi.formula = ~ DayAmdmt,
data = soil_phylum_small_otu1)
sand_vcov(mod)

}
