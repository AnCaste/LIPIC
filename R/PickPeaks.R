#' @title pick.peaks
#'
#' @description Function to identify local maxima in a vector, typically a spectrum or a
#' chromatogram.
#' N.B. This function belongs to the ChemometricsWithR package,
#' developed by Ron Wehrens.
#'
#' @param x Numerical vector.
#' @param span Neighbourhood, used to define local maxima.
#' 
#' @importFrom stats embed
#'
#' @return A vector containing positions of local maxima in the input data.
#' 
#' @references Ron Wehrens: Chemometrics with R. Springer Verlag, Berlin Hei-delberg (2020).
#' DOI: 10.1007/978-3-642-17841-2.
#'
#' @export

pick.peaks <- function(x, span) {
  span.width <- span * 2 + 1
  loc.max <- span.width + 1 -
    apply(embed(x, span.width), 1, which.max)
  loc.max[loc.max == 1 | loc.max == span.width] <- NA

  pks <- loc.max + 0:(length(loc.max)-1)
  unique(pks[!is.na(pks)])
}
