#' @title VisualSpec
#'
#' @description Allows the experimental mass spectra interactive visualization using the
#' "plotly" graphical interface. Moreover, spectral peak-picking is performed. VisualSpec
#' can be used to determine monoisotopic peak coordinates (m/z and Intensity) for a specific
#' lipid sum compositions.
#'
#' @param path File path leading to the Tab Separated Value (TSV) text file (.txt)
#' containing the m/z and Intensity values of each point of the experimental mass spectrum.
#' @param mzmin The lower boundary of the spectrum m/z range. mzmin=0 is set as default.
#' @param mzmax The upper boundary of the spectrum m/z range. mzmax=0 is set as default.
#' Specifically, when mzmax=0, the full experimental m/z range is visualized. Conversely,
#' if the user is interested in a specific spectral m/z range, the mzmin and mzmax
#' need to be specified.
#' @param span The "span" parameter value required by the "pick.peaks" function.
#' This function is included in the "LIPIC" package, but belongs to the
#' "ChemometricsWithR" package.
#'
#' @return A data frame containing the m/z and Intensity values related to peaks detected
#' in the mass spectrum by the "pick.peaks" function.
#' 
#' @import plotly
#' @importFrom utils read.delim
#'
#' @examples
#' \dontrun{
#'
#' > VisualSpec("C:/Documents/Spectra/CardiolipinMix.txt")
#' # The whole m/z range is considered.
#'   The file path ("C:/Documents/Spectra/CardiolipinMix.txt") is illustrative only.
#'   Please, replace it with your spectrum file (.txt) path.
#'
#' > VisualSpec("C:/Documents/Spectra/CardiolipinMix.txt",650,750)
#' # The 720-730 m/z spectral subrange is considered.
#'
#' }
#' 
#' @references Ron Wehrens: Chemometrics with R. Springer Verlag, Berlin Hei-delberg (2020).
#' DOI: 10.1007/978-3-642-17841-2.
#' 
#' @references Carson Sievert: "Interactive Web-Based Data Visualization with R, plotly, and shiny"
#' Chapman and Hall/CRC, 2020, https://plotly-r.com.
#'
#' @export

VisualSpec=function(path,mzmin=0,mzmax=0,span=25){
  spec=read.delim(path)
  if (mzmax==0){
    subspec=spec
  }else{
    subspec=spec[mzmin<spec[, 1] & spec[, 1]<mzmax,]
  }
  Spectrum=plot_ly(data=subspec,x=subspec[,1],y=subspec[,2],type='scatter',mode='lines',name="Spectrum",line = list(color = "rgba(0,0,255,0.4)"))
  Spectrum=Spectrum%>%layout(yaxis=list(title="Relative Abundance(%)",type='linear'),xaxis=list(title="m/z"))
  print(Spectrum)
  pn=pick.peaks(subspec[[2]],span)
  pmz=0
  pI=0
  for (i in 1:length(pn)){
    pmz[i]=subspec[pn[i],1]
    pI[i]=subspec[pn[i],2]
  }
  ODF=data.frame(pmz,pI)
  names(ODF)=c("Peak m/z", "Peak Intensity")
  return(ODF)
}
