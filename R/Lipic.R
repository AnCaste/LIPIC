#' @title Lipic
#'
#' @description Analogously to the "Lipic_It" function, "Lipic" performs both Type I and Type II
#' corrections in lipid class-related mass spectra. Furthermore, both "Lipic" and "Lipic_It" share
#' the same input and output (see the "Lipic_It" help page for further details).
#' However, the "Lipic" function is not based on a recursive workflow. Hence, the user should
#' refer to it only when reliable dmz values are somehow known.
#'
#' @param data A data frame reporting the recognized sum compositions in the first column.
#' Conversely, the corresponding chemical formulas, the measured m/z, the dmz, and the
#' absolute intensity values need to be listed in the second, third, fourth and fifth
#' columns respectively.
#' N.B. All the elements in the data frame must be sorted following the ascending
#' order of the measured m/z value. Moreover, the "dmz" column must be filled with the
#' corresponding estimates.
#' @param z Charge of the lipid ions.The charge sign needs to be specified.
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
#' The first five columns are identical to those of the input data frame.
#' The sixth column shows the intensity values after Type II correction. The seventh column shows
#' the intensity values after both Type II and Type I corrections. As indicated by the corresponding
#' headers, the remaining columns show the "Adj" and "Un-Adj" "Lipid Class" percent Relative
#' Abundance and the "Adj" and "Un-Adj" "Most Abundant Lipid" percent Relative Abundance. See the
#' description in the Lipic_It help page for details about the meaning of "Adj" and "Un-Adj".
#' 
#' @import plotly
#' @importFrom enviPat isopattern
#' @importFrom utils read.delim
#'
#' @examples
#' \dontrun{
#'
#' > DF=Lipic(data,-2,75000,span=25,ppm=10,xseq=0.001)
#' > Do you wish to compare the real spectrum with the simulated ones? (yes/no): yes
#' # Type yes if you need to check the Lipic correction reliability by
#'   simulated and real mass spectra comparison.Here, the "Lipic_It" output data frame
#'   is assigned to the "DF" variable.
#'
#' > Please, paste .txt mass spectrum file path: C:/Documents/Spectra/CLMixSpectrum.txt
#' # Paste the path leading to the .txt file related to the experimental spectrum.
#'   Do not write file path in quotes. Here, the file path
#'   ("C:/Documents/Spectra/CLMixSpectrum.txt") is illustrative only.
#'   Please, replace it with the file path leading to the position in which you want the
#'   .txt file to be saved.
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

Lipic=function(data,z,Res,t=0.001,span=25,ppm=10,xseq=0.001){
  
  isotopes=LIPIC::Isotope_Table

######## Lipic Core Algorithm ########

  names(data)=c("Sum Composition","Formula","Measured m/z","dmz","Intensity")
  n=dim(data)[1]
  IsoPattern=isopattern(isotopes,data[["Formula"]],threshold=t,charge=z,plotit="FALSE")
  sc=data[["Sum Composition"]]
  mz=data[["Measured m/z"]]
  dmz=data[["dmz"]]
  Ipc=data[["Intensity"]]
  Ic=0
  Inc=0
  Ic[1]=Ipc[1]*(sum(IsoPattern[[1]][,2])/IsoPattern[[1]][1,2])
  Inc[1]=Ic[1]
  for (i in 2:n){
    if (i==2){
      M1=IntFind(IsoPattern[[i-1]],IsoPattern[[i]][1,1],z)
      AL=sum(IsoPattern[[i]][,2])/IsoPattern[[i]][1,2]
      if (dim(M1)[1]>0){
        M1[,2]=M1[,2]/100
        c1=0
        for (j in 1:dim(M1)[1]){
          s=(M1[j,1]+dmz[i])/(Res*(2*sqrt(2*log(2))))
          c=M1[j,2]*exp(-(mz[i]-(M1[j,1]+dmz[i]))^2/(2*(s)^2))
          c1=c1+c
        }
        c1=c1*Ipc[i-1]
        s=(IsoPattern[[i]][1,1]+dmz[i])/(Res*(2*sqrt(2*log(2))))
        cL=exp(-(mz[i]-(IsoPattern[[i]][1,1]+dmz[i]))^2/(2*(s)^2))
        Inc[i]=Ipc[i]*AL
        Ipc[i]=((Ipc[i]-c1)/cL)
        if (Ipc[i]<0){
          Ic[i]=0
        }else{
          Ic[i]=Ipc[i]*AL
        }
      }else{
        Ic[i]=Ipc[i]*AL
        Inc[i]=Ic[i]
      }
    }else{
      M1=IntFind(IsoPattern[[i-1]],IsoPattern[[i]][1,1],z)
      M2=IntFind(IsoPattern[[i-2]],IsoPattern[[i]][1,1],z)
      AL=sum(IsoPattern[[i]][,2])/IsoPattern[[i]][1,2]
      if (dim(M1)[1]>0&dim(M2)[1]>0){
        M1[,2]=M1[,2]/100
        c1=0
        for (j in 1:dim(M1)[1]){
          s=(M1[j,1]+dmz[i])/(Res*(2*sqrt(2*log(2))))
          c=M1[j,2]*exp(-(mz[i]-(M1[j,1]+dmz[i]))^2/(2*(s)^2))
          c1=c1+c
        }
        M2[,2]=M2[,2]/100
        c2=0
        for (j in 1:dim(M2)[1]){
          s=(M2[j,1]+dmz[i])/(Res*(2*sqrt(2*log(2))))
          c=M2[j,2]*exp(-(mz[i]-(M2[j,1]+dmz[i]))^2/(2*(s)^2))
          c2=c2+c
        }
        c1=c1*Ipc[i-1]
        c2=c2*Ipc[i-2]
        s=(IsoPattern[[i]][1,1]+dmz[i])/(Res*(2*sqrt(2*log(2))))
        cL=exp(-(mz[i]-(IsoPattern[[i]][1,1]+dmz[i]))^2/(2*(s)^2))
        Inc[i]=Ipc[i]*AL
        Ipc[i]=((Ipc[i]-c1-c2)/cL)
        if (Ipc[i]<0){
          Ic[i]=0
        }else{
          Ic[i]=Ipc[i]*AL
        }
      }else if(dim(M1)[1]>0){
        M1[,2]=M1[,2]/100
        c1=0
        for (j in 1:dim(M1)[1]){
          s=(M1[j,1]+dmz[i])/(Res*(2*sqrt(2*log(2))))
          c=M1[j,2]*exp(-(mz[i]-(M1[j,1]+dmz[i]))^2/(2*(s)^2))
          c1=c1+c
        }
        c1=c1*Ipc[i-1]
        s=(IsoPattern[[i]][1,1]+dmz[i])/(Res*(2*sqrt(2*log(2))))
        cL=exp(-(mz[i]-(IsoPattern[[i]][1,1]+dmz[i]))^2/(2*(s)^2))
        Inc[i]=Ipc[i]*AL
        Ipc[i]=((Ipc[i]-c1)/cL)
        if (Ipc[i]<0){
          Ic[i]=0
        }else{
          Ic[i]=Ipc[i]*AL
        }
      }else{
        Ic[i]=Ipc[i]*AL
        Inc[i]=Ic[i]
      }
    }
  }

  PITc=(Ic/sum(Ic))*100
  ARc=(Ic/max(Ic))*100
  PITnc=(Inc/sum(Inc))*100
  ARnc=(Inc/max(Inc))*100
  updf=cbind(data,Ipc,Ic,PITc,PITnc,ARc,ARnc)
  names(updf)[5:11]=c("Measured Intensity","Type II","Type II & Type I","Lipid Class Adj-RA (%)","Lipid Class UnAdj-RA (%)","Most Abundant Lipid Adj-RA (%)","Most Abundant Lipid UnAdj-RA (%)")

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
    Spectra=Spectra%>%plotly::layout(yaxis=list(title="Relative Abundance(%)",type='linear'),xaxis=list(title="m/z"))
    print(Spectra)
  }
  return(updf)
}
