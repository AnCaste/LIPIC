#' @title dmzFind
#'
#' @description Provides an estimate of m/z peak shifts (dmz) in the real mass spectrum in
#' respect to the corresponding m/z values in the simulated (theoretical) spectrum obtained
#' after Type I and Type II corrections.
#'
#' @param sspec A data frame containing the coordinates
#' of points in the simulated mass spectrum.
#' @param dat A typical Lipic or Lipic_It input data frame.
#' @param span The "span" parameter value required by the "pick.peaks" function.
#' This function is included in the "LIPIC" package, but belongs to the "ChemometricsWithR"
#' package.
#' @param ppm Accuracy threshold (in ppm) within the m/z of a peak in the simulated mass
#' spectrum matches a measured m/z value among those listed in the third column of the
#' input "dat" data frame.
#'
#' @return A modified version of the input "dat" data frame in which the initial dmz values
#' are replaced by those calculated as the difference between the measured and the
#' theoretical m/z values which matched under the "ppm" threshold.
#' 
#' 
#' @references Ron Wehrens: Chemometrics with R. Springer Verlag, Berlin Hei-delberg (2020).
#' DOI: 10.1007/978-3-642-17841-2.

dmzFind=function(sspec,dat,span,ppm){
  x=pick.peaks(sspec[[2]],span)
  sp=0
  for (i in 1:length(x)){
    sp[i]=sspec[x[i],1]
  }
  r=dim.data.frame(dat)[1]
  mzm=dat[[3]]
  smzmax=replicate(r,0)
  for (j in 1:r){
    for (i in 1:length(sp)){
      if ((sp[i]-mzm[j])<((sp[i]*ppm)/1E6)){
        smzmax[j]=sp[i]
      }
    }
  }
  dmz=replicate(r,0)
  for (j in 1:r){
    if (smzmax[j]!=0){
      dmz[j]=mzm[j]-smzmax[j]
    }
  }
  dat[[4]]=dmz
  return(dat)
}
