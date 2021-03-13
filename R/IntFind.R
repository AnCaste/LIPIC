#' @title IntFind
#'
#' @description Detects isotopologues responsible for Type II isotopic effect.
#'
#' @param IM A matrix containing the m/z ratios (first column) and the Relative Abundances
#' (second column) of the interfering lipid isotopologues. These abundances are expressed
#' as percentage in respect to the preeminent isotopologue.
#' @param mzr m/z ratio of the M+0 isotopologue undergoing isotopic interference.
#' @param ch Charge of the ions involved in Type II isotopic effect.
#'
#' @return A matrix containing m/z and relative abundances of isotopologues causing the
#' Type II isotopic effect.

IntFind=function(IM,mzr,ch){
  m=0
  for (k in 1:length(IM[,1])){
    d=abs(mzr-(IM[k,1]))
    if (d<(0.4/abs(ch))){
      mv=c(IM[k,1],IM[k,2])
      m=c(m,mv)
    }
  }
  m=m[-1]
  m=matrix(m,ncol=2,byrow="TRUE")
  return(m)
}
