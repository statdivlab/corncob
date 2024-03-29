% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_taxa_names.R
\name{clean_taxa_names}
\alias{clean_taxa_names}
\title{Rename taxa}
\usage{
clean_taxa_names(x, name = "OTU")
}
\arguments{
\item{x}{Object of class \code{phyloseq}}

\item{name}{Character, defaults to \code{"OTU"}. Optional. String to use in every taxa name.}
}
\value{
Object of class \code{phyloseq}, with taxa renamed (defaults to OTU1, OTU2, ...), with the original taxa names saved as an attribute.
}
\description{
Renames taxa to have short human-readable names
}
\details{
The original taxa names are saved as the \code{original_names} attribute. See the example for an example of how to access the original names.
}
