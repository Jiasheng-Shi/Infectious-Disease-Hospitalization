library(mvtnorm)
library(EpiEstim)
library(mc2d)
library(ggplot2)
library(nlme)
library(numDeriv)
library(TTR)
library(parallel)
library(dplyr)
  
#####################################
# Real data, load precleaned data

{
  
  ### load the data
  #setwd()
  reg.dat<- read.csv('./derived_data.csv')
  #Miami-Dade 12086
  #Wayne MI 26163
  #NYC NY 36061
  #COOK IL 17031
  
}

#####################################
# Summarize of the reg.dat
{
  
  ### Summarize data set
  reg.dat.summarize <- reg.dat%>% group_by(fips)%>%
    summarize(county = unique(county), state= unique(state), 
              max_cases = max(inc,na.rm = TRUE),min_cases = min(inc,na.rm = TRUE),
              ob_days= sum(inc>0),
              min_ob_day=min(date),max_ob_day=max(date))
  View(reg.dat.summarize)
  na_flag <- apply(is.na(reg.dat), 2, sum)
  na_flag
  
}

#####################################
## Determine the length of two generation intervals
## eta1 and eta2
## Cross validation
## Together with seeking for initial values 

{
  
  names(reg.dat)[5]="H"
  
  ###########################################################################
  
  ### Composite for omega
  ### Fixed I and H
  
  V3compOmega<-function(paraT,omegafi,eta1,reg.dat){
    Rfi<-list()
    Lambdafi<-list()
    Td<-c()
    NumI=list()
    DenI=list()
    result=c()
    omegafi=c(omegafi[1:(eta1-1)],(1-sum(omegafi[1:(eta1-1)])))
    
    if (omegafi[eta1]<0){
      penalty=-omegafi[eta1]
      omegafi[eta1]=0
    }else{
      penalty=0
    }
    
    for (countyi in 1:nrow(reg.dat.summarize)) {
      ## Lambdafi
      reg.dat.inter=reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]
      Isave=reg.dat.inter$inc
      
      reg.dat.inter=reg.dat.inter[reg.dat.inter$date>=as.Date('2021-01-01','%Y-%m-%d') ,]
      
      Td[countyi]=nrow(reg.dat.inter)
      
      Lambdafi[[countyi]]<-WMA(c(Isave[(32-eta1):31],reg.dat.inter$inc),
                    n=eta1,wts = rev(omegafi))[eta1:(eta1+Td[countyi]-1)] 
      
      ## Rfi(R-function-intermediate)
      R0inter=Isave[32]/Isave[31]
      Rfi[[countyi]]=c(1:Td[countyi])
      Rfi[[countyi]][1]=paraT[1] + (paraT[2]*log(R0inter)) + as.numeric((as.numeric(reg.dat.inter[1,6:7])%*%paraT[3:4]))
      for (r in 2:Td[countyi]) {
        Rfi[[countyi]][r]= (paraT[1]) + (paraT[2]*log( reg.dat.inter$inc[(r-1)]/Lambdafi[[countyi]][(r-1)] ))+ 
          (as.numeric(reg.dat.inter[r,6:7])%*%paraT[3:4] )  
      }
      Rfi[[countyi]]=exp(Rfi[[countyi]])
      
      ## NumI & DenI
      NumI[[countyi]]=rep(0,Td[[countyi]])
      DenI[[countyi]]=rep(1,Td[[countyi]])
      for (i in 1:Td[[countyi]]) {
        NumI[[countyi]][i]=reg.dat.inter$inc[i]-Lambdafi[[countyi]][i]*Rfi[[countyi]][i]
        DenI[[countyi]][i]=Lambdafi[[countyi]][i]*Rfi[[countyi]][i]
      }
      
      NumI[[countyi]]=NumI[[countyi]]^2
      result[countyi]=sum(NumI[[countyi]]/DenI[[countyi]])
      
    }
    
    resultA=sum(result)+10000000*penalty
    
    return(resultA)
  }
  #
  
  V3compOmega2<-function(x){
    paraT=x[1:4]
    omegafi=x[5:(3+eta1)]
    result=V3compOmega(paraT,omegafi,eta1 = eta1,reg.dat=reg.dat)
    return(result)
  }
  
  simuResIni=list()
  
  
  timeS=Sys.time()
  for (eta1 in 22:22) {
    simuResIni[[eta1]]=optim(fn=V3compOmega2,par = c(0.6,0.5,-0.01,-0.12,rep(0.03,(eta1-1))), 
                     lower = c(c(-1,-1,-1,-1),rep(0.001,(eta1-1))),upper = c(c(1,1,1,1), rep(1,(eta1-1)) ),
                     method = "L-BFGS-B",
                     #control = list(pgtol=0 )
    ) 
  }
  Sys.time()-timeS
  
  ##Potential eta1=7,18,19,20,22,24
  
  ###########################################################################
  
  ### Composite for omega
  ### Fixed I and H
  
     reg.dat.TiH=reg.dat
     reg.dat.TiH=reg.dat.TiH%>%mutate(inc2 = lag(inc,n=1,default=0),
                              inc3 = lag(inc,n=2,default=0),
                              inc4 = lag(inc,n=3,default=0),
                              inc5 = lag(inc,n=4,default=0),
                              inc6 = lag(inc,n=5,default=0),
                              inc7 = lag(inc,n=6,default=0),
                              inc8 = lag(inc,n=7,default=0),
                              inc9 = lag(inc,n=8,default=0),
                              inc10 = lag(inc,n=9,default=0),
                              inc11 = lag(inc,n=10,default=0),
                              inc12 = lag(inc,n=11,default=0),
                              inc13 = lag(inc,n=12,default=0),
                              inc14 = lag(inc,n=13,default=0),
                              inc15 = lag(inc,n=14,default=0),
                              inc16 = lag(inc,n=15,default=0),
                              inc17 = lag(inc,n=16,default=0),
                              inc18 = lag(inc,n=17,default=0),
                              inc19 = lag(inc,n=18,default=0),
                              inc20 = lag(inc,n=19,default=0),
                              inc21 = lag(inc,n=20,default=0))
      
     reg.dat.TiH=reg.dat.TiH[reg.dat.TiH$date>=as.Date('2021-01-01','%Y-%m-%d') ,]
     
     V3compTildeOmega<-function(RETi,tildeomegafi,eta2,reg.dat.TiH){
       
       
       penalty=matrix(data = NA, nrow = nrow(reg.dat.summarize),ncol = eta2)
       reg.dat.TiH$Num=reg.dat.TiH$H
       
       for (countyi in 1:(nrow(reg.dat.summarize))) {
         if (RETi[countyi]<0){
           penalty[countyi,1]=-RETi[countyi]
           RETi[countyi]=0
         }else if(RETi[countyi]>1){
           penalty[countyi,1]=RETi[countyi]-1
           RETi[countyi]=1
         }else{
           penalty[countyi,1]=0
         }
         
         reg.dat.TiH[reg.dat.TiH$fips==reg.dat.summarize$fips[countyi],]$Num=
           reg.dat.TiH[reg.dat.TiH$fips==reg.dat.summarize$fips[countyi],]$Num-
           RETi[countyi]*reg.dat.TiH[reg.dat.TiH$fips==reg.dat.summarize$fips[countyi],][,8]
       }
         
       for (i in 1:(eta2-2)) {
         if (tildeomegafi[i]<0){
           penalty[,(i+1)]=rep((-tildeomegafi[i]),nrow(reg.dat.summarize))
           tildeomegafi[i]=0
         }else if(tildeomegafi[i]>1){
           penalty[,(i+1)]=rep((tildeomegafi[i]-1),nrow(reg.dat.summarize))
           tildeomegafi[i]=1
         }else{
           penalty[,(i+1)]=rep(0,nrow(reg.dat.summarize))
         }
           
         reg.dat.TiH$Num=reg.dat.TiH$Num-tildeomegafi[i]*reg.dat.TiH[,i+7]
           
       }
       
       for (countyi in 1:nrow(reg.dat.summarize)) {
         tildeomegafi[(eta2-1)]=1-RETi[countyi]-sum(tildeomegafi[1:(eta2-2)])
         if (tildeomegafi[(eta2-1)]<0){
           penalty[,eta2]=rep((-tildeomegafi[(eta2-1)]),nrow(reg.dat.summarize))
           tildeomegafi[(eta2-1)]=0
         }else if(tildeomegafi[(eta2-1)]>1){
           penalty[,eta2]=rep((tildeomegafi[(eta2-1)]-1),nrow(reg.dat.summarize))
           tildeomegafi[(eta2-1)]=1
         }else{
           penalty[,eta2]=rep(0,nrow(reg.dat.summarize))
         }
       }
         
       
         
         ## NumI & DenI
         reg.dat.TiH$Num=(reg.dat.TiH$Num)^2
         result=sum(reg.dat.TiH$Num)
       
       resultA=result+1e10*sum(penalty)
       
       return(resultA)
     }
     
     V3compTildeOmega2<-function(TransTildeo){
       RETi=rep(0,nrow(reg.dat.summarize))
       tildeomegafi=rep(0,(eta2-2))
       for (i in 1:(eta2-2)) {
         tildeomegafi[(i+1)]=sum(TransTildeo[i:(eta2-2)])
       }
       for (i in 1:nrow(reg.dat.summarize)) {
         RETi[i]=TransTildeo[(eta2-2+i)]+tildeomegafi[2]
       }
       result=V3compTildeOmega(RETi,tildeomegafi,eta2 = eta2,reg.dat.TiH=reg.dat.TiH)
       return(result)
     }
     #
      
     TilsimuResIni=list()
     
     timeS=Sys.time()
     for (eta2 in 5:5) {
       #TilsimuResIni[[eta2]]=optim(fn=V3compTildeOmega2,par = rep(0.045,(eta2-1)), 
       #                          lower = (abs(rnorm(n=(eta2-1),mean = 0.03, sd=0.01))),
        #                        upper =  rep(1,(eta2-1)),
        #                        method = "L-BFGS-B",
                                #control = list(pgtol=0 )
       #) 
       TilsimuResIni[[eta2]]=optim(fn=V3compTildeOmega2,par = rep(0.045,(eta2+2)), 
                                   lower = (abs(rnorm(n=(eta2+2),mean = 0.005, sd=0.001))),
                                   upper =  rep(1,(eta2+2)),
                                   method = "L-BFGS-B",
                                   #control = list(pgtol=0 )
       ) 
     }
     Sys.time()-timeS
  
  #
  
     OmegaIni=c(simuResIni[[22]]$par[5:25],1-sum(simuResIni[[22]]$par[5:25]))
     paraIni=simuResIni[[22]]$par[1:4]
     Tilinter=TilsimuResIni[[5]]$par
     tildeomegaIni=matrix(data = NA,nrow = nrow(reg.dat.summarize),ncol=5)
     for (i in 2:4) {
       tildeomegaIni[,i]=rep( sum(Tilinter[(i-1):3]),nrow(reg.dat.summarize))
     }
     for (countyi in 1:nrow(reg.dat.summarize)) {
       tildeomegaIni[countyi,1]=Tilinter[(3+countyi)]+tildeomegaIni[countyi,2]
       tildeomegaIni[countyi,5]=1-sum(tildeomegaIni[countyi,1:4])
     }
     eta1=22
     eta2=5
  
}

##################################
## Likelihood function with known generation intervals

{

  ## tuning parameter for how many sets of h to sample
  TN0=1000
  
  ## Function been invoked in Complike
  CompLFrbinom<-function(x){
    return(rbinom(n=TN0,size = ceiling(x[1]),prob = x[2]))
  }
  
  CompLFdbinom<-function(xx){
    return(dbinom(x=xx[1],size = xx[2],prob = xx[3],log = TRUE))
  }
  
  ########### Main Composite Likelihood function
  CompLike<-function(paraTOutput,paraT,omegaOutput,tildeomegaOutput,omegaInput,tildeomegaInput){
    ## Rfi(R-function-intermediate) and Lambdafi
    Rfi<-list()
    Lambdafi<-list()
    RfiOut<-c()
    LambdafiOut<-c()
    Td<-c()
    #hfi=list()
    #cellQlog=list()
    #TimeQlog=c()
    
    reg.dat$Lambdafi=NA
    reg.dat$Rfi=NA
    reg.dat$LambdafiOut=NA
    reg.dat$RfiOut=NA
    reg.dat$LogLike=NA
    
    for (countyi in 1:nrow(reg.dat.summarize)) {
      
      hfi=list()
      cellQlog=list()
      
      ## Lambdafi
      
      reg.dat.inter=reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]
      Isave=reg.dat.inter$inc
      
      reg.dat.inter=reg.dat.inter[reg.dat.inter$date>=as.Date('2021-01-01','%Y-%m-%d') ,]
      Hsave=reg.dat.inter$H
      
      Td[countyi]=nrow(reg.dat.inter)
      
      Lambdafi[[countyi]]<-WMA(c(Isave[(32-eta1):31],reg.dat.inter$inc),
                               n=eta1,wts = rev(omegaInput))[eta1:(eta1+Td[countyi]-1)] 
      
      ## Rfi(R-function-intermediate)
      R0inter=Isave[32]/Isave[31]
      Rfi[[countyi]]=c(1:Td[countyi])
      Rfi[[countyi]][1]=paraT[1] + (paraT[2]*log(R0inter)) + as.numeric((as.numeric(reg.dat.inter[1,6:7])%*%paraT[3:4]))
      for (r in 2:Td[countyi]) {
        Rfi[[countyi]][r]= (paraT[1]) + (paraT[2]*log( reg.dat.inter$inc[(r-1)]/Lambdafi[[countyi]][(r-1)] ))+ 
          (as.numeric(reg.dat.inter[r,6:7])%*%paraT[3:4] )  
      }
      Rfi[[countyi]]=exp(Rfi[[countyi]])
      
      reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]$Lambdafi[32:(31+Td[countyi])]=Lambdafi[[countyi]]
      reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]$Rfi[32:(31+Td[countyi])]=Rfi[[countyi]]

      ##############################################
      
      Iforh=c(Isave[(34-eta2):31],Isave[32:(31+Td[countyi])])
      H=reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]$H[32:(31+Td[countyi])]
      for (r in 1:Td[countyi]) {
        pforEh=cbind(Iforh[(r+length(tildeomegaInput[countyi,])-3):r],tildeomegaInput[countyi,2:(eta2-1)])
        pforEh<-lapply(seq_len(nrow(pforEh)), function(i) pforEh[i,])
        
        hfi[[r]]=mapply(FUN = CompLFrbinom, pforEh)
        
        hfi[[r]]=cbind( ( rep(Hsave[r],TN0)-rowSums(hfi[[r]]) ) ,hfi[[r]])
        hfi[[r]]=hfi[[r]][hfi[[r]][,1]>=0,]
        
        successP=dpois(x=hfi[[r]][,1],lambda = (tildeomegaInput[countyi,1]*Lambdafi[[countyi]][r]*Rfi[[countyi]][r]),log = TRUE)
        
        successP=successP+dpois(x=max(0,(ceiling(Isave[(r+31)])-hfi[[r]][,1])),lambda = ((1-tildeomegaInput[countyi,1])*Lambdafi[[countyi]][r]*Rfi[[countyi]][r]),log = TRUE )
        
        successP=log(2)+successP+log(pi)+ log(Lambdafi[[countyi]][r])+log(Rfi[[countyi]][r])+ log( tildeomegaInput[countyi,1]*(1-tildeomegaInput[countyi,1]) )/2
        
        successP=successP- mean(successP)
        successP[successP<=-30]=-Inf
        while (mean(successP[successP>-30])>5) {
          successP=successP-mean(successP[successP>=-30])
          successP[successP<=-30]=-Inf
        }
        
        successP=exp(successP)
        successP=successP/(max(successP) )
        successP[successP>1]=1
        
        successPsi=mapply(FUN = rbinom, n=1,size=1, prob=successP)
        
        hfi[[r]]=cbind(hfi[[r]],successPsi)
        
        if (sum(successPsi)<=(TN0/100)){
          rowindex=which(hfi[[r]][,length(tildeomegaInput[countyi,])]==1)
          hfi[[r]]=hfi[[r]][c(1:(TN0/200), rep(rowindex,(TN0/100))),]
        }else{
          hfi[[r]]=hfi[[r]][hfi[[r]][,length(tildeomegaInput[countyi,])]==1,] 
        }
        
        hfi[[r]]=hfi[[r]][,-length(tildeomegaInput[countyi,])]
      }
      
      LambdafiOut[[countyi]]<-WMA(c(Isave[(32-eta1):31],reg.dat.inter$inc),
                       n=eta1,wts = rev(omegaOutput))[eta1:(eta1+Td[countyi]-1)]
      
      RfiOut[[countyi]]=c(1:Td[countyi])
      RfiOut[[countyi]][1]= paraTOutput[1] + (paraTOutput[2]*log(R0inter)) + as.numeric((as.numeric(reg.dat.inter[1,6:7])%*%paraTOutput[3:4]))
      for (r in 2:Td[countyi]) {
        RfiOut[[countyi]][r]= paraTOutput[1] + (paraTOutput[2]*log( reg.dat.inter$inc[(r-1)]/Lambdafi[[countyi]][(r-1)] ))+
          (as.numeric(reg.dat.inter[r,6:7])%*%paraTOutput[3:4] ) 
      }
      RfiOut[[countyi]]=exp(RfiOut[[countyi]])
      
      reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]$LambdafiOut[32:(31+Td[countyi])]=LambdafiOut[[countyi]]
      reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]$RfiOut[32:(31+Td[countyi])]=RfiOut[[countyi]]
      
      CompLFdbinomdpois<-function(xx,rr){
        y=c()
        y[1]=dpois(x=xx[1],lambda=(tildeomegaOutput[countyi,1]*RfiOut[[countyi]][rr]*LambdafiOut[[countyi]][rr]),log = TRUE)
        y[2]=dpois(x=max(0,(ceiling(Isave[(rr+31)])-xx[1])),lambda = ((1-tildeomegaOutput[countyi,1])*RfiOut[[countyi]][rr]*LambdafiOut[[countyi]][rr]),log = TRUE)
        yy=cbind(xx[-1],ceiling(Iforh[(rr+eta2-3):rr]),tildeomegaOutput[countyi,2:(eta2-1)])
        yy<-lapply(seq_len(nrow(yy)), function(i) yy[i,])
        y[3:eta2]=mapply(FUN = CompLFdbinom, yy)
        return(y)
      }
      
      ## Calculate the \hat{Q}
      for (r in 1:Td[countyi]) {
        
        interhForQ=lapply(seq_len(nrow(hfi[[r]])), function(i) hfi[[r]][i,])
        
        cellQlog[[r]]=t(mapply(FUN=CompLFdbinomdpois, interhForQ ,rr=r))
      }
      
      cellQoverall=mapply(sum, cellQlog)
      cellQN=mapply(nrow, cellQlog)
      reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]$LogLike[32:(31+Td[countyi])]=cellQoverall/cellQN
        
    }
    
    cellQLV=sum(na.omit(reg.dat$LogLike))
    
    resCompLike=list(cellQLV=cellQLV,reg.dat=reg.dat)
    return(resCompLike)
    
  }
  
}

########## Negative CompLike function with fixed omega and tildeomega
########## First Iteration, paraIni initial values

{
    
    NCLFOmega1=function(paraTOutput){
      resNCL=-CompLike(paraTOutput,paraT=paraIni,
                       omegaOutput=OmegaIni,
                       tildeomegaOutput=tildeomegaIni,
                       omegaInput=OmegaIni,
                       tildeomegaInput=tildeomegaIni)$cellQLV
      #resNCL=-sum(resNCL$cellQoverall/resNCL$cellQN)
      return(resNCL)
    }  
  
  paraIni
  
  timeS=Sys.time()
  simuRes1=nlm(NCLFOmega1,p= paraIni,steptol = 1e-5)
  Sys.time()-timeS
  
  simuRes1
  
  V3compOmega3<-function(x){
    omegafi=x[1:(eta1-1)]
    result=V3compOmega(paraT=simuRes1$estimate,omegafi,eta1 = eta1,reg.dat=reg.dat)
    return(result)
  }
  
  OmegaIni
  
  timeS=Sys.time()
  simuomega1=optim(fn=V3compOmega3,par = OmegaIni[1:21], 
                   lower = rep(0.001,21),upper = rep(1,21),
                   method = "L-BFGS-B",
                   #control = list(pgtol=0 )
  )
  Sys.time()-timeS
  
  simuomega1$par
  Omega1=c(simuomega1$par,1-sum(simuomega1$par))
  paraT1=simuRes1$estimate
  
  tildeomega1=tildeomegaIni

}

########## Negative CompLike function with fixed omega and tildeomega
########## Second Iteration, arbitrary initial values

{

  NCLFOmega2=function(paraTOutput){
    resNCL=-CompLike(paraTOutput,paraT=simuRes1$estimate,
                     omegaOutput=Omega1,
                     tildeomegaOutput=tildeomega1,
                     omegaInput=Omega1,
                     tildeomegaInput=tildeomega1)$cellQLV
    #resNCL=-sum(resNCL$cellQoverall/resNCL$cellQN)
    return(resNCL)
  }
  
  paraT1
  
  timeS=Sys.time()
  simuRes2=nlm(NCLFOmega2,p= simuRes1$estimate,steptol = 1e-6)
  Sys.time()-timeS
  
  simuRes2
  
  V3compOmega4<-function(x){
    omegafi=x[1:(eta1-1)]
    result=V3compOmega(paraT=simuRes2$estimate,omegafi,eta1 = eta1,reg.dat=reg.dat)
    return(result)
  }
  
  Omega1
  
  timeS=Sys.time()
  simuomega2=optim(fn=V3compOmega4,par = Omega1[1:21], 
                   lower = rep(0.001,21),upper = rep(1,21),
                   method = "L-BFGS-B",
                   #control = list(pgtol=0 )
  )
  Sys.time()-timeS
  
  simuomega2$par
  Omega2=c(simuomega2$par,1-sum(simuomega2$par))
  paraT2=simuRes2$estimate
  
  tildeomega2=tildeomegaIni

}

##################################

#saveRDS(list(simuRes1,simuRes2,simuRes3,simuResIni,I,H,Z,simuomega1,simuomega2,simuomega3,tildeomegaIni), file='HosCovidCase3.RDS')
  
  RealDataRes=list(simuRes1,simuRes2,simuResIni,reg.dat,TilsimuResIni,tildeomegaIni,eta1,eta2)
  

#Sys.time()-timerep

saveRDS(RealDataRes,file="RealDataResult.RDS")


##################################
# Bootstrap for omega

###########################################################################

### Composite for omega
### Fixed I and H

timeS=Sys.time()
{

  Bitrep=200
  BootOmega=list()
  BootSample=list()
  for (Boot in 1:(Bitrep)) {
    BootSample[[Boot]]=sample(181,90,replace = FALSE)
    
    BootV3compOmega<-function(paraT,omegafi,eta1,reg.dat){
      Rfi<-list()
      Lambdafi<-list()
      Td<-c()
      NumI=list()
      DenI=list()
      result=c()
      omegafi=c(omegafi[1:(eta1-1)],(1-sum(omegafi[1:(eta1-1)])))
      
      if (omegafi[eta1]<0){
        penalty=-omegafi[eta1]
        omegafi[eta1]=0
      }else{
        penalty=0
      }
      
      for (countyi in 1:nrow(reg.dat.summarize)) {
        ## Lambdafi
        reg.dat.inter=reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]
        Isave=reg.dat.inter$inc
        
        reg.dat.inter=reg.dat.inter[reg.dat.inter$date>=as.Date('2021-01-01','%Y-%m-%d') ,]
        
        Td[countyi]=nrow(reg.dat.inter)
        
        Lambdafi[[countyi]]<-WMA(c(Isave[(32-eta1):31],reg.dat.inter$inc),
                                 n=eta1,wts = rev(omegafi))[eta1:(eta1+Td[countyi]-1)] 
        
        ## Rfi(R-function-intermediate)
        R0inter=Isave[32]/Isave[31]
        Rfi[[countyi]]=c(1:Td[countyi])
        Rfi[[countyi]][1]=paraT[1] + (paraT[2]*log(R0inter)) + as.numeric((as.numeric(reg.dat.inter[1,6:7])%*%paraT[3:4]))
        for (r in 2:Td[countyi]) {
          Rfi[[countyi]][r]= (paraT[1]) + (paraT[2]*log( reg.dat.inter$inc[(r-1)]/Lambdafi[[countyi]][(r-1)] ))+ 
            (as.numeric(reg.dat.inter[r,6:7])%*%paraT[3:4] )  
        }
        Rfi[[countyi]]=exp(Rfi[[countyi]])
        
        ## NumI & DenI
        NumI[[countyi]]=rep(0,Td[[countyi]])
        DenI[[countyi]]=rep(1,Td[[countyi]])
        for (i in 1:Td[[countyi]]) {
          NumI[[countyi]][i]=reg.dat.inter$inc[i]-Lambdafi[[countyi]][i]*Rfi[[countyi]][i]
          DenI[[countyi]][i]=Lambdafi[[countyi]][i]*Rfi[[countyi]][i]
        }
        
        NumI[[countyi]]=NumI[[countyi]]^2
        result[countyi]=sum( (NumI[[countyi]]/DenI[[countyi]])[BootSample[[Boot]]] )
        
      }
      
      resultA=sum(result)+10000000*penalty
      
      return(resultA)
    }
    #
    
    BootV3compOmega2<-function(x){
      omegafi=x[1:(eta1-1)]
      result=BootV3compOmega(paraT=paraT2,omegafi,eta1 = eta1,reg.dat=reg.dat)
      return(result)
    }
    
      BootOmega[[Boot]]=optim(fn=BootV3compOmega2,par = c(rep(0.03,(eta1-1))), 
                               lower = rep(0.001,(eta1-1)),upper = rep(1,(eta1-1)),
                               method = "L-BFGS-B")
    
  }
  
  for (boot in 1:Bitrep) {
    
    BootOmega[[boot]]$par=c(BootOmega[[boot]]$par,1-sum(BootOmega[[boot]]$par))
    
  }
  
}

Sys.time()-timeS

BootResSave=list(BootSample,BootOmega)
saveRDS(BootResSave,file="BootResSave.RDS")

#########################################################################
############################################ Omega Voronoi Treemap

{

  library(voronoiTreemap)
  library(latex2exp)
  library(RColorBrewer)
  library(cowplot)
  
  data(ExampleGDP)
  vt_input_from_df(ExampleGDP)
  gdp_json <- vt_export_json(vt_input_from_df(ExampleGDP))
  
  ExampleGDP$h2="Asia"
  ExampleGDP=ExampleGDP[1:22,]
  ExampleGDP$weight=Omega2
  #ExampleGDP$color[1:22] <- scales::seq_gradient_pal(high = "#928087",low = "#F4F3EB")(ExampleGDP$weight[1:22]/max(ExampleGDP$weight[1:22]))
  
  ExampleGDP$color[1:22] <- scales::seq_gradient_pal(high = "#895A58",low = "#F4F3EB")(ExampleGDP$weight[1:22]/max(ExampleGDP$weight[1:22]))
  gdp_json <- vt_export_json(vt_input_from_df(ExampleGDP,scaleToPerc = FALSE))
  vt_d3(gdp_json,size_circle = "3px", color_circle = "white",legend = FALSE,seed = 1100,label = FALSE)
  
  #652*409
}

#########################################################################
############################################  TildeOmega Voronoi Treemap

{

  data(ExampleGDP)
  vt_input_from_df(ExampleGDP)
  gdp_json <- vt_export_json(vt_input_from_df(ExampleGDP))
  
  ExampleGDP$h2="Asia"
  ExampleGDP=ExampleGDP[1:5,]
  ExampleGDP$weight=tildeomega2[4,]
  ExampleGDP$color[1:5]<-c(brewer.pal(6,"PuRd")[5:2],"#ece9ea")
  gdp_json <- vt_export_json(vt_input_from_df(ExampleGDP,scaleToPerc = FALSE))
  vt_d3(gdp_json,size_circle = "3px", color_circle = "white",legend = FALSE,seed = 10,label = FALSE)
  
  #588*417
  
}

#########################################################################
############################################
# Bootstrap CI for four county estimated R

simuResBootPara=list()
timeS=Sys.time()
for (boot in 1:Bitrep) {
  
  NCLFBoot=function(paraTOutput){
    resNCL=-CompLike(paraTOutput,paraT=paraIni,
                     omegaOutput=BootOmega[[boot]]$par,
                     tildeomegaOutput=tildeomegaIni,
                     omegaInput=BootOmega[[boot]]$par,
                     tildeomegaInput=tildeomegaIni)$cellQLV
    #resNCL=-sum(resNCL$cellQoverall/resNCL$cellQN)
    return(resNCL)
  }  
  
  simuResBootPara[[boot]]=nlm(NCLFBoot,p= paraIni,steptol = 1e-5)
  
}
Sys.time()-timeS

BootResSave=list(BootSample,BootOmega,simuResBootPara)
saveRDS(BootResSave,file="BootResSave.RDS")


#########################################################################
###Four county R 

{
  
  EstR<-list()
  reg.dat$EstR<-NA
  Td=c()
  for (countyi in 1:nrow(reg.dat.summarize)) {
    
    reg.dat.inter=reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]
    Isave=reg.dat.inter$inc
    
    reg.dat.inter=reg.dat.inter[reg.dat.inter$date>=as.Date('2021-01-01','%Y-%m-%d') ,]
    Hsave=reg.dat.inter$H
    
    Td[countyi]=nrow(reg.dat.inter)
    
    ## Rfi(R-function-intermediate)
    R0inter=Isave[32]/Isave[31]
    EstR[[countyi]]=c(1:Td[countyi])
    EstR[[countyi]][1]=paraT2[1] + (paraT2[2]*log(R0inter)) + as.numeric((as.numeric(reg.dat.inter[1,6:7])%*%paraT2[3:4]))
    for (r in 2:Td[countyi]) {
      EstR[[countyi]][r]= (paraT2[1]) + (paraT2[2]* EstR[[countyi]][(r-1)] )+ 
        (as.numeric(reg.dat.inter[r,6:7])%*%paraT2[3:4] )  
    }
    EstR[[countyi]]=exp(EstR[[countyi]])
    
    reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]$EstR[32:(31+Td[countyi])]=EstR[[countyi]]
    
    ##############################################
    
  }
  
  BootEstR=list()
  for (countyi in 1:nrow(reg.dat.summarize)) {
    
    BootEstR[[countyi]]=matrix(data = NA, nrow=Td[countyi], ncol=10)
    
    reg.dat.inter=reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]
    Isave=reg.dat.inter$inc
    
    reg.dat.inter=reg.dat.inter[reg.dat.inter$date>=as.Date('2021-01-01','%Y-%m-%d') ,]
    Hsave=reg.dat.inter$H
    
    R0inter=Isave[32]/Isave[31]
    
    for (boot in 1:Bitrep) {
      
      interP=simuResBootPara[[boot]]$estimate
      ## Rfi(R-function-intermediate)
      BootEstR[[countyi]][1,boot]=interP[1]+interP[2]*log(R0inter)+
        as.numeric((as.numeric(reg.dat.inter[1,6:7])%*%interP[3:4]))
      for (r in 2:Td[countyi]) {
        BootEstR[[countyi]][r,boot]= (interP[1]) + (interP[2]* BootEstR[[countyi]][(r-1),boot] )+ 
          (as.numeric(reg.dat.inter[r,6:7])%*%interP[3:4] )  
      }
      BootEstR[[countyi]][,boot]=exp(BootEstR[[countyi]][,boot])
      
    }
    
  }
  
  reg.dat$BL1=NA
  reg.dat$BL2=NA
  reg.dat$BU1=NA
  reg.dat$BU2=NA
  for (countyi in 1:nrow(reg.dat.summarize)) {
    for (r in 1:Td[countyi]) {
      reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]$BL1[(r+31)]=min(BootEstR[[countyi]][r,])
      reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]$BU1[(r+31)]=max(BootEstR[[countyi]][r,])
      reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]$BL2[(r+31)]=sort(BootEstR[[countyi]][r,])[2]
      reg.dat[reg.dat$fips==reg.dat.summarize$fips[countyi],]$BU2[(r+31)]=sort(BootEstR[[countyi]][r,])[9]
    }
  }
  
  reg.dat.gg=reg.dat[reg.dat$date>=as.Date('2021-01-01','%Y-%m-%d') ,]
  
  reg.dat.gg$BL1=reg.dat.gg$BL1-rep(runif(4,min=0,max=0.05),each=nrow(reg.dat.gg)/4)
  reg.dat.gg$BU2=reg.dat.gg$BU2+rep(runif(4,min=0,max=0.05),each=nrow(reg.dat.gg)/4)
  
  reg.dat.gg <-rbind(reg.dat.gg[,c(2,3,9,10,13)], reg.dat.gg[,c(2,3,9,10,13)])
  reg.dat.gg[(1+nrow(reg.dat.gg)/2):(nrow(reg.dat.gg)), 3] <- rep(1,(nrow(reg.dat.gg)/2))
  reg.dat.gg$Line <- rep(c('Estimated Instantaneous Reproduction Number', 'Critical Value R=1'), each=(nrow(reg.dat.gg)/2))
  reg.dat.gg$date=rep(c(1:181),8)
  
  reg.dat.gg[reg.dat.gg$county=="Cook",]$county="Cook, IL"
  reg.dat.gg[reg.dat.gg$county=="Miami-Dade",]$county="Miami-Dade, FL"
  reg.dat.gg[reg.dat.gg$county=="Wayne",]$county="Wayne, MI"
  reg.dat.gg[reg.dat.gg$county=="New York",]$county="New York, NY"
  
  
  
  yuplim=1.16
  ylowlim=0.75
  ggplot(reg.dat.gg, aes(date)) + 
    geom_line(aes(y=EstR,colour=Line, lty=Line)) + 
    #geom_line(aes(y=1), colour="black",lty=2)+
    labs(
      x = "Date",
      y = "instantaneous reproduction number"
    )+
    #geom_ribbon(aes(ymin=BL1, ymax=BU2),fill="#88a0b8",alpha=0.41)+
    #geom_ribbon(aes(ymin=BL1, ymax=BU2),fill="#d2aec9",alpha=0.41)+
    geom_ribbon(aes(ymin=BL1, ymax=BU2),fill="gray",alpha=0.53)+
    facet_grid(county ~ .)+
    theme_bw() + 
    theme(legend.position="top")+
    # theme(
    #   panel.grid.minor = element_blank(),
    #   legend.title = element_blank(),
    #   legend.position = c(0.75,0.9) 
    # ) +
    annotate(
      geom = "rect",
      ymin = -Inf, ymax = Inf, 
      xmin = 91,
      xmax = 121,
      #fill = c("#b8b8d0"), alpha = 0.6)+
      #fill = c("#08ffc4"), alpha = 0.15)+
      fill = c("#99cc00"), alpha = 0.45)+
    scale_color_manual(values=c("dodgerblue", "#CC0033"))+
    #scale_color_manual(values=c("#405880","orchid"))+
    #scale_color_manual(values=c("dodgerblue", "black"))+
    scale_linetype_manual(values=c("twodash", "solid")) +
    #ylim(ylowlim, (yuplim)) +
    scale_x_continuous(name='Date', breaks=c(1, 32, 60, 91, 121, 152,182), 
                       labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul"), limits=c(0,182) ) 
  
  #save, 1045*660

}

#########################################################################
### Omega Shape

library(ggridges)
{
  
  OmegaBoot=matrix(data = NA, nrow = eta1, ncol = 20)
  for (boot in 1:20) {
    OmegaBoot[,boot]=BootOmega[[boot]]$par
  }
  
  OmegaPoint=data.frame(TSI=c(1:22),Omega=Omega2)
  OmegaPoint$Up=NA
  OmegaPoint$Low=NA
  for (i in 1:eta1) {
    OmegaPoint$Up[i]=sort(OmegaBoot[i,])[19]
    OmegaPoint$Low[i]=sort(OmegaBoot[i,])[2]
  }
  
  OmegaPic=c(0,Omega2)
  OmegaCurve=spline(x=c(0:22),y=OmegaPic,n=5*length(OmegaPic))
  OmegaCurve=data.frame(TSI=OmegaCurve$x,Omega=OmegaCurve$y)
  OmegaCurve$Up=spline(x=c(0:22),y=c(0,OmegaPoint$Up),n=5*length(OmegaPic))$y
  OmegaCurve$Low=spline(x=c(0:22),y=c(0,OmegaPoint$Low),n=5*length(OmegaPic))$y
  
  Arrow1 <- data.frame(x1 = 10.7, x2 = 8.7, y1 = 0.07, y2 = 0.06)
  
  ggplot(OmegaCurve,aes(x=TSI))+
    geom_line(aes(x=TSI,y=Omega),lty=4)+
    geom_ribbon_pattern(
      aes(ymin=0, ymax=Omega), 
      pattern = 'image', 
      pattern_filename = "virus.png",
      pattern_type     = 'tile',
      pattern_scale    = -2,
    #  fill    = 'white'
    ) +
    #scale_pattern_fill_gradient(
    #  low = "red",
    #  high = "blue",
      #space = "Lab",
      #na.value = "grey50",
      #guide = guide_colourbar(available_aes = "pattern_fill"),
      #aesthetics = "pattern_fill"
    #)+
    #scale_pattern_type_manual(values = c("hexagonal", "rhombille","pythagorean"))+
    #geom_ribbon(data=OmegaCurve,aes(ymin=0, ymax=Omega) , fill="#d0c8d0",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve,aes(ymin=0, ymax=Omega) , fill="steelblue",alpha=0.25)+
    geom_ribbon(data=OmegaCurve,aes(ymin=0, ymax=Omega) , fill="white",alpha=0.6)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=2,],aes(ymin=0, ymax=Omega) , fill="#b06870",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=4,],aes(ymin=0, ymax=Omega) , fill="#b87880",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=6,],aes(ymin=0, ymax=Omega) , fill="#e09090",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=8,],aes(ymin=0, ymax=Omega) , fill="#d098b0",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=10,],aes(ymin=0, ymax=Omega) , fill="#f0b0b8",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=12,],aes(ymin=0, ymax=Omega) , fill="#e0b0b8",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=14,],aes(ymin=0, ymax=Omega) , fill="#d8a8b0",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=16,],aes(ymin=0, ymax=Omega) , fill="#f8d8d8",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI>=15.7,],aes(ymin=0, ymax=Omega) , fill="#f0d0d8",alpha=0.5)+
    
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=2,],aes(ymin=0, ymax=Omega) , fill="#f8d8d8",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=4 & OmegaCurve$TSI>=1.7,],aes(ymin=0, ymax=Omega) , fill="#f0d0d8",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=8 & OmegaCurve$TSI>=5.7,],aes(ymin=0, ymax=Omega) , fill="#f0b0b8",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=10 & OmegaCurve$TSI>=7.7,],aes(ymin=0, ymax=Omega) , fill="#e0b0b8",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=12 & OmegaCurve$TSI>=9.7,],aes(ymin=0, ymax=Omega) , fill="#b06870",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=14 & OmegaCurve$TSI>=11.7,],aes(ymin=0, ymax=Omega) , fill="#b87880",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI<=16 & OmegaCurve$TSI>=13.7,],aes(ymin=0, ymax=Omega) , fill="#d8a8b0",alpha=0.5)+
    #geom_ribbon(data=OmegaCurve[OmegaCurve$TSI>=15.7,],aes(ymin=0, ymax=Omega) , fill="#d098b0",alpha=0.5)+
    xlab("time since infection (day)")+
    ylab("infectiousness (density)")+
    geom_point(data=OmegaPoint,mapping = aes(x = TSI, y = Omega),
               colour="#ea3c53",size=2)+
    #theme_bw()+
    #theme(
       #panel.grid.minor = element_blank(),
     #  legend.title = element_blank(),
      # legend.position = c(15,0.09) 
    #)+
    geom_curve(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      data = Arrow1,
      curvature = 0.6,
      arrow = arrow(length = unit(0.03, "npc"))
    )+
    annotate(geom = "text",
             x = 15.3, y = 0.07,
             label = "infectiousness profile (viral shedding)",
             size=3.75)
    
  # 652*417
  
}

#########################################################################
### Tildeomega Shape

{
  
  TOmegaPoint=data.frame(TSI=rep(c(1:4),4),TOmega=c(t(tildeomega2[,1:4])))
  TOmegaPoint$group=rep(c("Miami-Dade, FL","Cook, IL", "Wayne, MI", "New York, NY"),each=4)
  
  #OmegaPic=c(0,Omega2)
  TOmegaCurve1=spline(x=c(0,c(1:4)-0.5,4.5),y=c(0,TOmegaPoint$TOmega[1:4],0),n=5*length(OmegaPic))
  TOmegaCurve2=spline(x=c(0,c(1:4)-0.5,4.5),y=c(0,TOmegaPoint$TOmega[5:8],0),n=5*length(OmegaPic))
  TOmegaCurve3=spline(x=c(0,c(1:4)-0.5,4.5),y=c(0,TOmegaPoint$TOmega[9:12],0),n=5*length(OmegaPic))
  TOmegaCurve4=spline(x=c(0,c(1:4)-0.5,4.5),y=c(0,TOmegaPoint$TOmega[13:16],0),n=5*length(OmegaPic))
  TOmegaCurve=data.frame(x=c(TOmegaCurve1$x,TOmegaCurve2$x,TOmegaCurve3$x,TOmegaCurve4$x),
                         y=c(TOmegaCurve1$y,TOmegaCurve2$y,TOmegaCurve3$y,TOmegaCurve4$y),
                         group=rep(c("Miami-Dade, FL","Cook, IL", "Wayne, MI", "New York, NY"),each=115) )
  
  TOmegaPoint$TSI=TOmegaPoint$TSI-0.5
  
  ggplot(TOmegaPoint)+
    geom_col(aes(x=TSI,y=TOmega,col=group,alpha=0.5))+
    geom_line(data = TOmegaCurve,aes(x=x,y=y,col=group))+
    facet_grid(~ group)+
    theme(legend.position = "none")+
    xlab("Time (day)")+
    ylab("hospitalization profile")
    
  # 652*417
  
}







