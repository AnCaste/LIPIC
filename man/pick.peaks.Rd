% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PickPeaks.R
\name{pick.peaks}
\alias{pick.peaks}
\title{pick.peaks}
\usage{
pick.peaks(x, span)
}
\arguments{
\item{x}{Numerical vector.}

\item{span}{Neighbourhood, used to define local maxima.}
}
\value{
A vector containing positions of local maxima in the input data.
}
\description{
Function to identify local maxima in a vector, typically a spectrum or a
chromatogram.
N.B. This function belongs to the ChemometricsWithR package,
developed by Ron Wehrens.
}
\references{
Ron Wehrens: Chemometrics with R. Springer Verlag, Berlin Hei-delberg (2020).
DOI: https://doi.org/10.1007/978-3-642-17841-2.
}
