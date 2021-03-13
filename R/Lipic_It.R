#' @title Lipic_It
#'
#' @description Performs both Type I and Type II corrections in lipid class-related mass spectra.
#'  To this purpose, however, the dmz values (i.e. the absolute error in monoisotopic m/z
#'  estimation) need to be specified for each sum composition.
#'  Typically, such values are initially unknown. In this case, the user should employ
#'  "Lipic_It" rather than "Lipic" function. Indeed, "Lipic_It" follows a recursive algorithm
#'  leading to a reliable estimation of dmz values. After the last iteration, the Type I and Type II
#'  correction results are shown.
#'  Specifically, "Lipic_It" returns the values of the lipid Relative Abundances compared with both
#'  the "Most Abundant Lipid" and the "Lipid Class".
#'  The former indicates the percentage ratio between the signal intensity of a specific
#'  lipid sum composition and the highest intensity detected in the mass spectrum.
#'  The latter corresponds to the percentage ratio between the signal intensity of
#'  that sum composition and the overall intensity measured for the whole lipid class.
#'  If equal ionization yields are assessed, "Lipid Class" percent Relative Abundances
#'  can be considered equivalent to the mole fraction percentage of the sum compositions.
#'  In the "Lipic_It" output data frame, both "Most Abundant Lipid" and "Lipid Class"
#'  percent Relative Abundances are defined as "Adj" or "Un-Adj". Precisely, "Adj"
#'  (Adjusted) indicates that the relative abundances are computed after both Type I
#'  and Type II corrections. Conversely, Un-Adj (UnAdjusted) Relative Abundances are
#'  calculated using only Type I correction.
#'  The bar chart generated in the RStudio viewer pane at the end of "Lipic_It"
#'  calculations shows the interactive comparison between the "Adj" and "Un-Adj"
#'  values of "Lipid Class" percent Relative Abundances. This graphical output allows the
#'  user to perceive the extent of Type II correction in his experimental mass spectrum.
#'  Moreover, the reliability of the corrections performed by "Lipic_It" can be evaluated
#'  by requesting the comparison between the real mass spectrum and the
#'  simulated mass spectra obtained after either Type I or both Type I and Type II
#'  corrections. Subsequently, a new plot is generated in the RStudio viewer pane.
#'  Here, these spectra are named "Uncorrected Simulated Spectrum" (Type I correction) and
#'  "Corrected Simulated Spectrum" (Type I and Type II correction). The user can also
#'  choose which spectra need to be removed from spectral superimposition by clicking on
#'  the corresponding line in the legend panel.
#'
#' @param data A data frame reporting the recognized sum compositions in the first column.
#' Conversely, the corresponding chemical formulas, the measured m/z, the dmz, and the
#' absolute intensity values need to be listed in the second, third, fourth and fifth
#' columns respectively.
#' N.B. All the elements in the data frame must be sorted following the ascending
#' order of the measured m/z value. Furthermore, the "dmz" column must be filled with 0.
#' Please, check the "CLmix" data frames included in the LIPIC package
#' as an example of data organization of the "Lipic_It" input data frame.
#' @param z Charge of the lipid ions. The charge sign needs to be specified.
#' @param Res Operating resolution of the mass analyzer.
#' @param t Probability below which isotope peaks can be omitted.
#' Refer to the "isopattern" function in the "enviPat" package for details.
#' @param span The "span" parameter value required by the "pick.peaks" function.
#' This function is included in the "LIPIC" package, but belongs to the
#' "ChemometricsWithR" package.
#' @param ppm Accuracy threshold (in ppm) required by the "dmzFind" function.
#' @param xseq Increment in the sequence of the m/z values used for the calculation
#' of the simulated mass spectrum employed in dmz calculation. "xseq" is set to 0.001
#' as default. Higher values reduce the computation time, but can lower the accuracy by
#' which the dmz values are retrieved by the "dmzFind" function.
#'
#' @return A data frame which is an extension of the input data frame.
#' The first five columns are identical to those of the input data frame, except the dmz
#' column (i.e. the fourth column), which contains the dmz estimated values for each sum composition.
#' The accuracy of "Lipic_It" correction is strictly dependent on the reliability of dmz values.
#' These are recursively calculated until convergence and the overall content of the output data
#' frame reflects the results of the last iteration.
#' Precisely, the sixth column shows the intensity values after Type II correction.
#' The seventh column shows the intensity values after both TypeII and Type I corrections.
#' As indicated by the corresponding headers, the remaining columns show the "Adj" and "Un-Adj"
#' "Lipid Class" percent Relative Abundance and the "Adj" and "Un-Adj" "Most Abundant Lipid"
#' percent Relative Abundance (see the above description section for details).
#' 
#' @import plotly
#' @importFrom utils read.delim
#'
#' @examples
#' \dontrun{
#'
#' # Follow these instructions to perform the "Lipic_It" isotope correction on MS data
#' referred to a standard Cardiolipin mix purchased from Avanti Polar Lipids.
#'
#' > data(CLmix)
#' # Load the "CLmix" input data frame.
#'
#' > data(CLMixSpectrum)
#' # Load the experimental mass spectrum of the cardiolipin mix.
#'
#' > write.table(CLMixSpectrum,"C:/Documents/Spectra/CLMixSpectrum.txt",sep="\t",row.names=FALSE)
#' # Export the "CLMixSpectrum" data frame in a Tab Separated Value text file (.txt).
#'   The file path ("C:/Documents/Spectra/CLMixSpectrum.txt") is illustrative only.
#'   Please, replace it with the file path leading to the position in which you want the
#'   .txt file to be saved.
#'
#' > ODF=Lipic_It(CLmix,-2,75000,span=25,ppm=10,xseq=0.001)
#' > Do you wish to compare the real spectrum with the simulated ones? (yes/no): yes
#' # As shown, type yes if you need to check the Lipic_It correction reliability by
#'   comparison of Corrected and Uncorrected simulated spectra comparison with
#'   the real mass spectrum. Here, the "Lipic_It" output data frame is assigned to the
#'   "ODF" variable.
#'
#' > Please, paste .txt mass spectrum file path: C:/Documents/Spectra/CLMixSpectrum.txt
#' # Paste the path leading to the .txt file related to the cardiolipin mix mass
#'   spectrum. As shown in the example above, do not write the file path in quotes.
#'
#' }
#'
#' @references Andrea Castellaneta, Ilario Losito, Davide Coniglio, Beniamino Leoni,
#' Pietro Santamaria, Maria A. Di Noia, Luigi Palmieri, Cosima D. Calvano, Tommaso R.I. Cataldi:
#' "LIPIC: an automated workflow to account for isotopologues-related interferences in
#' electrospray ionization high resolution mass spectra of phospholipids",
#' Journal of the American Society for Mass Spectrometry, 2021. DOI: https://doi.org/10.1021/jasms.1c00008.
#' 
#' @references Martin Loos, Christian Gerber, Francesco Corona, Juliane Hollender,
#' Heinz Singer: "Accelerated Isotope Fine Structure Calculation Using Pruned Transition Trees"
#' Anal. Chem. 87, 5738–5744 (2015). DOI: https://doi.org/10.1021/acs.analchem.5b00941.
#' 
#' @references Ron Wehrens: Chemometrics with R. Springer Verlag, Berlin Hei-delberg (2020).
#' DOI: 10.1007/978-3-642-17841-2.
#' 
#' @references Carson Sievert: "Interactive Web-Based Data Visualization with R, plotly, and shiny"
#' Chapman and Hall/CRC, 2020, https://plotly-r.com.
#'
#' @export

Lipic_It=function(data,z,Res,t=0.001,span=25,ppm=10,xseq=0.001){
  n=dim(data)[1]
  lDF=iLipic(data,z,Res,t,span,ppm,xseq)
  DF=lDF[[1]]
  data[[4]]=DF[[4]]
  DF2=iLipic(data,z,Res,t,span,ppm,xseq)
  DF2=DF2[[1]]
  data[[4]]=DF2[[4]]
  c=0
  for (i in 1:n){
    if (DF[i,4]==DF2[i,4]){
      c=c+1
    }
  }
  while(c!=n){
    lDF=iLipic(data,z,Res,t,span,ppm,xseq)
    DF=lDF[[1]]
    data[[4]]=DF[[4]]
    DF2=iLipic(data,z,Res,t,span,ppm,xseq)
    DF2=DF2[[1]]
    data[[4]]=DF2[[4]]
    c=0
    for (i in 1:n){
      if (DF[i,4]==DF2[i,4]){
        c=c+1
      }
    }
  }

  ### Variable Assignment for Graphical Output ###

  IsoPattern=lDF[[2]]
  Ic=lDF[[3]]
  Inc=lDF[[4]]
  PITc=(Ic/sum(Ic))*100
  PITnc=(Inc/sum(Inc))*100
  sc=data[[1]]

  ######## Lipic Graphical Part: Bar Charts ########

  MI=as.character(PITnc)
  AI=as.character(PITc)
  Hist=plot_ly(y=MI, x=sc, type = "bar",name="Lipid class unadjusted relative Abundance (%)",marker = list(color = "rgba(0, 0, 255, 0.4)"))
  Hist=Hist %>% add_trace(y=AI, x=sc, type = "bar",name="Lipid class adjusted relative Abundance (%)",marker = list(color = "rgba(204,204,0, 0.7)"))
  Hist=Hist %>% plotly::layout(yaxis=list(title="Abundance(%)",type='linear'),xaxis=list(title="Sum composition",categoryorder = "array",categoryarray =sc),bargap=0.6)
  print(Hist)

  ######## Lipic Graphical Part: Spectral Overlap ########

  Q=readline(prompt="Do you wish to compare the real spectrum with the simulated ones? (yes/no): ")
  if (Q=="yes"){
    path=readline(prompt="Please, paste .txt mass spectrum file path: ")
    up=IsoPattern[[n]][length(IsoPattern[[n]][,1]),1]+2
    lo=IsoPattern[[1]][1,1]-2
    df=read.delim(path)
    names(df)=c("m/z","Relative Abundance (%)")
    sdf=df[lo<df[, 1] & df[, 1]< up,c("m/z","Relative Abundance (%)")]
    x=seq(lo,up,0.001)
    yc=0
    ync=0
    for (j in 1:n){
      s=(IsoPattern[[j]][,1])/(Res*(2*sqrt(2*log(2))))
      ys=0
      yn=0
      for (i in 1:length(IsoPattern[[j]][,1])){
        pc=((Ic[j]*IsoPattern[[j]][i,2])/sum(IsoPattern[[j]][,2]))*exp(-(x-IsoPattern[[j]][i,1])^2/(2*s[i]^2))
        pnc=((Inc[j]*IsoPattern[[j]][i,2])/sum(IsoPattern[[j]][,2]))*exp(-(x-IsoPattern[[j]][i,1])^2/(2*s[i]^2))
        ys=ys+pc
        yn=yn+pnc
      }
      yc=yc+ys
      ync=ync+yn
    }
    yc=(yc/max(sdf[,2]))*100
    ync=(ync/max(sdf[,2]))*100
    sdf[,2]=(sdf[,2]/max(sdf[,2]))*100
    Spectra=plot_ly(x =x,y=yc,type='scatter',mode='lines',name="Corrected Simulated Spectrum",line = list(color = "rgba(0,0,255,0.4)"))
    Spectra=Spectra%>%add_trace(x=x,y=ync,type='scatter',mode='lines',name="Uncorrected Simulated Spectrum",line = list(color = "rgba(255,0,0,0.5)"))
    Spectra=Spectra%>%add_trace(data=sdf,x=sdf[,1],y=sdf[,2],type='scatter',mode='lines',name="Real Spectrum",line = list(color = "rgba(204,204,0,0.7)"))
    Spectra=Spectra%>%layout(yaxis=list(title="Relative Abundance(%)",type='linear'),xaxis=list(title="m/z"))
    print(Spectra)
  }
  return(DF)
}
