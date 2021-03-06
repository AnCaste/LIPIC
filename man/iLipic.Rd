% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iLipic.R
\name{iLipic}
\alias{iLipic}
\title{iLipic}
\usage{
iLipic(data, z, Res, t = 0.001, span = 25, ppm = 10, xseq = 0.001)
}
\arguments{
\item{data}{This parameter is automatically set by the "Lipic_It" function.}

\item{z}{This parameter is automatically set by the "Lipic_It" function.}

\item{Res}{This parameter is automatically set by the "Lipic_It" function.}

\item{t}{This parameter is automatically set by the "Lipic_It" function.}

\item{span}{This parameter is automatically set by the "Lipic_It" function.}

\item{ppm}{This parameter is automatically set by the "Lipic_It" function.}

\item{xseq}{This parameter is automatically set by the "Lipic_It" function.}
}
\value{
A list of elements useful to build the Lipic_It output.
}
\description{
iLipic is a computational tool invoked by the "Lipic_It" function for the recursive
correction of both Type I and Type II isotopic effects.
}
\references{
Andrea Castellaneta, Ilario Losito, Davide Coniglio, Beniamino Leoni,
Pietro Santamaria, Maria A. Di Noia, Luigi Palmieri, Cosima D. Calvano, Tommaso R.I. Cataldi:
"LIPIC: an automated workflow to account for isotopologues-related interferences in
electrospray ionization high resolution mass spectra of phospholipids",
Journal of the American Society for Mass Spectrometry, 2021. DOI: https://doi.org/10.1021/jasms.1c00008.

Martin Loos, Christian Gerber, Francesco Corona, Juliane Hollender,
Heinz Singer: "Accelerated Isotope Fine Structure Calculation Using Pruned Transition Trees"
Anal. Chem. 87, 5738–5744 (2015). DOI: https://doi.org/10.1021/acs.analchem.5b00941.

Ron Wehrens: Chemometrics with R. Springer Verlag, Berlin Hei-delberg (2020).
DOI: https://doi.org/10.1007/978-3-642-17841-2.
}
