#' @title iLipic
#'
#' @description iLipic is a computational tool invoked by the "Lipic_It" function for the recursive
#' correction of both Type I and Type II isotopic effects.
#'
#' @param data This parameter is automatically set by the "Lipic_It" function.
#' @param z This parameter is automatically set by the "Lipic_It" function.
#' @param Res This parameter is automatically set by the "Lipic_It" function.
#' @param t This parameter is automatically set by the "Lipic_It" function.
#' @param span This parameter is automatically set by the "Lipic_It" function.
#' @param ppm This parameter is automatically set by the "Lipic_It" function.
#' @param xseq This parameter is automatically set by the "Lipic_It" function.
#'
#' @return A list of elements useful to build the Lipic_It output.
#' 
#' @importFrom enviPat isopattern
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

iLipic=function(data,z,Res,t=0.001,span=25,ppm=10,xseq=0.001){
  
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

  #Simulated spectrum calculation

  up=IsoPattern[[n]][length(IsoPattern[[n]][,1]),1]+2
  lo=IsoPattern[[1]][1,1]-2
  x=seq(lo,up,xseq)
  yc=0
  for (j in 1:n){
    s=(IsoPattern[[j]][,1])/(Res*(2*sqrt(2*log(2))))
    ys=0
    for (i in 1:length(IsoPattern[[j]][,1])){
      pc=((Ic[j]*IsoPattern[[j]][i,2])/sum(IsoPattern[[j]][,2]))*exp(-(x-IsoPattern[[j]][i,1])^2/(2*s[i]^2))
      ys=ys+pc
    }
    yc=yc+ys
  }
  yc=(yc/max(yc))*100
  sspec=data.frame(x,yc)
  data2=dmzFind(sspec,data,span,ppm)

  #Calculation on new dataset

  names(data2)=c("Sum Composition","Formula","Measured m/z","dmz","Intensity")
  n=dim(data2)[1]
  IsoPattern=isopattern(isotopes,data2[["Formula"]],threshold=t,charge=z,plotit="FALSE")
  sc=data2[["Sum Composition"]]
  mz=data2[["Measured m/z"]]
  dmz=data2[["dmz"]]
  Ipc=data2[["Intensity"]]
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
  updf=cbind(data2,Ipc,Ic,PITc,PITnc,ARc,ARnc)
  names(updf)[5:11]=c("Measured Intensity","Type II","Type II & Type I","Lipid Class Adj-RA (%)","Lipid Class UnAdj-RA (%)","Most Abundant Lipid Adj-RA (%)","Most Abundant Lipid UnAdj-RA (%)")

  l=list(updf,IsoPattern,Ic,Inc)
  return(l)
}
