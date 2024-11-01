library(mvtnorm)
library(EpiEstim)
library(mc2d)
#library(ggplot2)
library(nlme)
#library(numDeriv)
library(TTR)
#library(parallel)
library(dplyr)

args <- commandArgs(TRUE)  
irep <- as.numeric(args[length(args)])

#timerep=Sys.time()
OResRep=list()

for (Simrep in 1:100) {
    
  # Generating data for simulation purpose (Correct model specification)
  
  ###################################
  ## Generate generation intervals
  
  {
    
    ## tildeomega[k]=precentage of patients went to the hospital after (k-1)days 
    ##    after his infection, k-1=0,...5
    ## tildeomega[7]=precentage of patients never went to the hospital after his infection, 
    tildeomega=c()
    tildeomega[7]=0.5
    for (i in 1:5) {
      tildeomega[i]=pgamma(i,shape = 1.6,scale = 1.5)-pgamma((i-1),shape = 1.6,scale = 1.5)
    }
    tildeomega[6]=pgamma(5,shape = 1.6,scale = 1.5,lower.tail = FALSE)
    tildeomega[1:6]=tildeomega[1:6]*(1-tildeomega[7])
    
    plot(x=c(1:6),y=tildeomega[1:6],type="l",lty=4)
    
    ## omega[k]=precentage of patients get infected by a infectious one after k days, k=1,...25
    omega=c()
    for (i in 1:24) {
      omega[i]=pgamma(i,shape = 2.5,scale = 3)-pgamma((i-1),shape = 2.5,scale = 3)
    }
    omega[25]=pgamma(24,shape = 2.5,scale = 3,lower.tail = FALSE)
    
    plot(x=c(1:25),y=omega,type="l",lty=4)
    
    
  }
  
  ###################################
  ## Generate the instantaneous reproduction number 
  
  {

    ## length of the Time series/number of observations
    T=120
    ## Generate Covariants
    ## Z is the covariates, Z[,1:(NoCov=2)] assumed to be time variate
    NoCov=2
    
    Case4LoadData=readRDS(file = "Case4LoadData.RDS")
    
    Z=Case4LoadData[[1]]
    
    ##################################
    ## True parameter value:
    ## Oracle beta
    OracleBeta=as.matrix(c(-0.02,-0.125),nrow=NoCov,ncol=1)
    
    ## Oracle phi
    OraclePhi=c(0.7,0.5)
    OraclePhi_0=OraclePhi[1]
    OraclePhi_1=OraclePhi[2]
    
    ##################################
    ## Generate the instantaneous reproduction number 
    ## R[t], the instantaneous reproduction number at time t.
    R<-matrix(data = NA, nrow=1,ncol=T)
    ### R[0], the instantaneous reproduction number at time 0.
    R_0=2.5
    
    #epsilon<-rep(0,T)
    R[1]=exp( (OraclePhi_0) + (OraclePhi_1*log(R_0)) + (Z[1,]%*%OracleBeta)  )
    for (t in 2:T) {
      R[t]=exp( (OraclePhi_0) + (OraclePhi_1*log(R[(t-1)]))+ (Z[t,]%*%OracleBeta)  )
    }
    plot(R[,],type = "l",lty=4)
    lines(x=c(0:120),y=rep(1,121),type="l",lty=4)
    
  }
  
  ##################################
  ## Generate the oracle incident cases
  
  {
    
    I_0=200
    
    ## Generate I[t], number of incidence at time t.
    I<-matrix(data = NA, nrow=1,ncol=T)
    I[1]=rpois(1,lambda = (R[1]*I_0*omega[1]))
    for (t in 2:25) { I[t]=rpois(1,lambda = ( R[t]*((I_0*omega[t])+ (I[1:(t-1)]%*%omega[(t-1):1] )  )  ) ) }
    for (t in 26:T) {
      if (( I[(t-1):(t-25)]%*%omega[1:25]  )<=100000){
        I[t]=rpois(1,lambda = (R[t]*( I[(t-1):(t-25)]%*%omega[1:25]  ) ))
      }else {
        multipconstant= (R[t]*( I[(t-1):(t-25)]%*%omega[1:25]  ) )%/%100000  
        residueconstant= (R[t]*( I[(t-1):(t-25)]%*%omega[1:25]  ) )%%100000  
        I[t]=sum( rpois(multipconstant,lambda = 100000   ) ) + rpois(1,lambda = residueconstant  )
      }
    }
    
    Isave=I
    
  }
  
  ##################################
  ## Generate the oracle hospital admission cases
  
  {
    
    h=rmultinomial(n=T,size = I, prob = tildeomega)
    H=matrix(data=NA,nrow=1,ncol=T)
    for (i in 1:T) {
      H[i]=0
      for (j in max(1,i+2-length(tildeomega)):i ) {
        H[i]=H[i]+h[j,(i+1-j)]
      }
    }
    
  }
  
  ##################################
  ## Missingness in incident cases
  
  {
    
    Missingrate=rnorm(T, mean = 0.15, sd = 0.05)
    Missingrate=abs(Missingrate)
    # for (i in 1:T) {
    #   if (i%%7!=0 && i%%7!=6){
    #     Missingrate[i]=0
    #   }
    # }
    
    Imiss=floor(I*(1-Missingrate))
   
    I=Imiss
    Isave
    
  }
  
  ##################################
  ## Seeking for initial values
  {
    ###########################################################################
    
    ### Composite for omega
    ### Fixed I and H
    
    V3compOmega<-function(paraT,omegafi,I,Z){
      Rfi<-c()
      Lambdafi<-c()
      NumI=c()
      DenI=c()
      omegafi=c(omegafi,(1-sum(omegafi)))
      
      if (omegafi[25]<0){
        penalty=-omegafi[25]
        omegafi[25]=0
      }else{
        penalty=0
      }
      
      ## Rfi(R-function-intermediate) and Lambdafi
      Lambdafi<-WMA(c(rep(0,(length(omega)-1)),I_0,c(I)),
                    n=length(omega),
                    wts = rev(omegafi))[length(omega):(length(omega)+T-1)]
      
      Rfi[1]=paraT[1] + (paraT[2]*log(R_0)) + (Z[1,]%*%paraT[3:(2+ncol(Z))])
      for (r in 2:T) {
        Rfi[r]= (paraT[1]) + (paraT[2]*log( I[r-1]/Lambdafi[r-1] ))+ (Z[r,]%*%paraT[3:(2+ncol(Z))] )  
      }
      Rfi=exp(Rfi)
      
      ## NumI & DenI
      for (i in 1:length(omega)) {
        NumI[i]=0
        DenI[i]=1
      }
      for (i in (length(omega)+1):T) {
        NumI[i]=I[i]-Lambdafi[i]*Rfi[i]
        DenI[i]=Lambdafi[i]*Rfi[i]
      }
      
      NumI=NumI^2
      result=sum(NumI/DenI)
      
      result=result+10000000*penalty
      
      return(result)
    }
    
    V3compOmega2<-function(x){
      paraT=x[1:4]
      omegafi=x[5:28]
      result=V3compOmega(paraT,omegafi,I=I,Z=Z)
      return(result)
    }
    
    omega
    
    timeS=Sys.time()
    simuResIni=optim(fn=V3compOmega2,par = c(0.6,0.5,-0.01,-0.12,rep(0.03,24)), 
                     lower = c(c(-1,-1,-1,-1),rep(0.001,24)),upper = c(c(1,1,1,1), rep(1,24) ),
                     method = "L-BFGS-B",
                     #control = list(pgtol=0 )
    )
    Sys.time()-timeS
    
    simuResIni$par
    OmegaIni=c(simuResIni$par[5:28],1-sum(simuResIni$par[5:28]))
    paraIni=simuResIni$par[1:4]
    
    ###########################################################################
    
    ### Composite for omega
    ### Fixed I and H
    
    reg.dat=data.frame(H=t(H))
    reg.dat$I1=c(I)
    mutate(reg.dat)
    
    reg.dat=reg.dat%>%mutate(I2 = lag(I1,n=1,default=0), 
                             I3 = lag(I1,n=2,default=0),
                             I4 = lag(I1,n=3,default=0),
                             I5 = lag(I1,n=4,default=0),
                             I6 = lag(I1,n=5,default=0))
    
    reg.dat=reg.dat[6:nrow(reg.dat),]
    
    tildeomegaIni=lm(H~ 0+I1+I2+I3+I4+I5+I6 ,data=reg.dat, )
    tildeomegaIni=c(tildeomegaIni$coefficients,1-sum(tildeomegaIni$coefficients))
    
    
  }
  
  ##################################
  ## Likelihood function with known generation intervals
  
  {
    
    ## tuning parameter for how many sets of h to sample
    TN0=1000
    
    ## Function been invoked in Complike
    CompLFrbinom<-function(x){
      return(rbinom(n=TN0,size = x[1],prob = x[2]))
    }
    
    CompLFdbinom<-function(xx){
      return(dbinom(x=xx[1],size = xx[2],prob = xx[3],log = TRUE))
    }
    
    ########### Main Composite Likelihood function
    CompLike<-function(paraTOutput,paraT,I,H,Z,omegaOutput,tildeomegaOutput,omegaInput,tildeomegaInput){
      ## Rfi(R-function-intermediate) and Lambdafi
      Rfi<-c()
      Lambdafi<-c()
      RfiOut<-c()
      LambdafiOut<-c()
      hfi=list()
      cellQlog=list()
      TimeQlog=c()
      
      ## Rfi(R-function-intermediate) and Lambdafi
      Rfi[1]=paraT[1] + (paraT[2]*log(R_0)) + (Z[1,]%*%paraT[3:(2+ncol(Z))])
      for (r in 2:T) {
        Rfi[r]= (paraT[1]) + (paraT[2]*Rfi[(r-1)])+ (Z[r,]%*%paraT[3:(2+ncol(Z))] )  
      }
      Rfi=exp(Rfi)
      
      Lambdafi<-WMA(c(rep(0,(length(omegaInput)-1)),I_0,c(I)),
                    n=length(omegaInput),
                    wts = rev(omegaInput))[length(omegaInput):(length(omegaInput)+T-1)]
      
      ## Sample small h, TN0 is a tuning parameter for how many sets of h to sample
      
      Iforh=c(rep(0,length(tildeomegaInput)-2),I)
      for (r in 1:T) {
        pforEh=cbind(Iforh[(r+length(tildeomegaInput)-3):r],tildeomegaInput[2:(length(tildeomegaInput)-1)])
        pforEh<-lapply(seq_len(nrow(pforEh)), function(i) pforEh[i,])
        
        hfi[[r]]=mapply(FUN = CompLFrbinom, pforEh)
        
        hfi[[r]]=cbind( ( rep(H[r],TN0)-rowSums(hfi[[r]]) ) ,hfi[[r]])
        hfi[[r]]=hfi[[r]][hfi[[r]][,1]>=0,]
        
        successP=dpois(x=hfi[[r]][,1],lambda = (tildeomegaInput[1]*Lambdafi[r]*Rfi[r]),log = TRUE)
        
        successP=successP+dpois(x=(I[r]-hfi[[r]][,1]),lambda = ((1-tildeomegaInput[1])*Lambdafi[r]*Rfi[r]),log = TRUE )
        
        successP=log(2)+successP+log(pi)+ log(Lambdafi[r])+log(Rfi[r])+ log( tildeomegaInput[1]*(1-tildeomegaInput[1]) )/2
        
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
          rowindex=which(hfi[[r]][,length(tildeomegaInput)]==1)
          hfi[[r]]=hfi[[r]][c(1:(TN0/200), rep(rowindex,(TN0/100))),]
        }else{
          hfi[[r]]=hfi[[r]][hfi[[r]][,length(tildeomegaInput)]==1,] 
        }
        
        hfi[[r]]=hfi[[r]][,-length(tildeomegaInput)]
      }
      
      
      ## RfiOut(R-function-intermediate-for-output, i.e., the unknown parameter is using omegaOutput 
      ## and tildeomegaOutput and paraTOutput) and LambdafiOut
      RfiOut[1]= paraTOutput[1] + (paraTOutput[2]*log(R_0)) + (Z[1,]%*%paraTOutput[3:(2+ncol(Z))])
      for (r in 2:T) {
        RfiOut[r]= paraTOutput[1] + (paraTOutput[2]*RfiOut[(r-1)])+ (Z[r,]%*%paraTOutput[3:(2+ncol(Z))] ) 
      }
      RfiOut=exp(RfiOut)
      
      LambdafiOut<-WMA(c(rep(0,(length(omegaInput)-1)),I_0,c(I)),
                       n=length(omegaInput),
                       wts = rev(omegaOutput))[length(omegaInput):(length(omegaInput)+T-1)]
      
      CompLFdbinomdpois<-function(xx,rr){
        y=c()
        y[1]=dpois(x=xx[1],lambda=(tildeomegaOutput[1]*RfiOut[rr]*LambdafiOut[rr]),log = TRUE)
        y[2]=dpois(x=(I[rr]-xx[1]),lambda = ((1-tildeomegaOutput[1])*RfiOut[rr]*LambdafiOut[rr]),log = TRUE)
        yy=cbind(xx[-1],Iforh[(rr+length(tildeomegaInput)-3):rr],tildeomegaOutput[2:(length(tildeomegaInput)-1)])
        yy<-lapply(seq_len(nrow(yy)), function(i) yy[i,])
        y[3:length(tildeomegaInput)]=mapply(FUN = CompLFdbinom, yy)
        return(y)
      }
      
      
      ## Calculate the \hat{Q}
      for (r in 1:T) {
        
        interhForQ=lapply(seq_len(nrow(hfi[[r]])), function(i) hfi[[r]][i,])
        
        cellQlog[[r]]=t(mapply(FUN=CompLFdbinomdpois, interhForQ ,rr=r))
      }
      
      cellQoverall=mapply(sum, cellQlog)
      cellQN=mapply(nrow, cellQlog)
      weightcell=rep(1,T)
      weightcell[which(cellQN<50)]=0
      cellQLV=sum(cellQoverall*weightcell/cellQN)
      
      resCompLike=list(cellQLV=cellQLV,cellQoverall=cellQoverall,cellQN=cellQN,cellQlog=cellQlog)
      return(resCompLike)
      
    }
    
  }
  
  ########## Negative CompLike function with fixed omega and tildeomega
  ########## First Iteration, paraIni initial values
  
  {
    
    NCLFOmega1=function(paraTOutput){
      resNCL=-CompLike(paraTOutput,paraT=paraIni,
                       I=I,
                       H=H,
                       Z=Z,
                       omegaOutput=OmegaIni,
                       tildeomegaOutput=tildeomegaIni,
                       omegaInput=OmegaIni,
                       tildeomegaInput=tildeomegaIni)$cellQLV
      #resNCL=-sum(resNCL$cellQoverall/resNCL$cellQN)
      return(resNCL)
    }
    
    c(OraclePhi,OracleBeta)
    
    timeS=Sys.time()
    simuRes1=nlm(NCLFOmega1,p= paraIni,steptol = 1e-5)
    Sys.time()-timeS
    
    simuRes1
    
    V3compOmega3<-function(x){
      omegafi=x[1:24]
      result=V3compOmega(paraT=simuRes1$estimate,omegafi,I=I,Z=Z)
      return(result)
    }
    
    omega
    
    timeS=Sys.time()
    simuomega1=optim(fn=V3compOmega3,par = OmegaIni[1:24], 
                     lower = rep(0.001,24),upper = rep(1,24),
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
                       I=I,
                       H=H,
                       Z=Z,
                       omegaOutput=Omega1,
                       tildeomegaOutput=tildeomega1,
                       omegaInput=Omega1,
                       tildeomegaInput=tildeomega1)$cellQLV
      #resNCL=-sum(resNCL$cellQoverall/resNCL$cellQN)
      return(resNCL)
    }
    
    c(OraclePhi,OracleBeta)
    
    timeS=Sys.time()
    simuRes2=nlm(NCLFOmega2,p= simuRes1$estimate,steptol = 1e-6)
    Sys.time()-timeS
    
    V3compOmega4<-function(x){
      omegafi=x[1:24]
      result=V3compOmega(paraT=simuRes2$estimate,omegafi,I=I,Z=Z)
      return(result)
    }
    
    omega
    
    timeS=Sys.time()
    simuomega2=optim(fn=V3compOmega4,par = OmegaIni[1:24], 
                     lower = rep(0.001,24),upper = rep(1,24),
                     method = "L-BFGS-B",
                     #control = list(pgtol=0 )
    )
    Sys.time()-timeS
    
    simuomega2$par
    Omega2=c(simuomega2$par,1-sum(simuomega2$par))
    paraT2=simuRes2$estimate
    
    tildeomega2=tildeomegaIni
    
  }
  
  ########## Negative CompLike function with fixed omega and tildeomega
  ########## Third Iteration, arbitrary initial values
  
  {
    
    NCLFOmega3=function(paraTOutput){
      resNCL=-CompLike(paraTOutput,paraT=simuRes2$estimate,
                       I=I,
                       H=H,
                       Z=Z,
                       omegaOutput=Omega2,
                       tildeomegaOutput=tildeomega2,
                       omegaInput=Omega2,
                       tildeomegaInput=tildeomega2)$cellQLV
      #resNCL=-sum(resNCL$cellQoverall/resNCL$cellQN)
      return(resNCL)
    }
    
    c(OraclePhi,OracleBeta)
    
    timeS=Sys.time()
    simuRes3=nlm(NCLFOmega3,p= simuRes2$estimate,steptol = 1e-6)
    Sys.time()-timeS
    
    V3compOmega5<-function(x){
      omegafi=x[1:24]
      result=V3compOmega(paraT=simuRes3$estimate,omegafi,I=I,Z=Z)
      return(result)
    }
    
    omega
    
    timeS=Sys.time()
    simuomega3=optim(fn=V3compOmega5,par = OmegaIni[1:24], 
                     lower = rep(0.001,24),upper = rep(1,24),
                     method = "L-BFGS-B",
                     #control = list(pgtol=0 )
    )
    Sys.time()-timeS
    
    simuomega3$par
    Omega3=c(simuomega3$par,1-sum(simuomega3$par))
    paraT3=simuRes3$estimate
    
    tildeomega3=tildeomegaIni
    
  }
  
  ##################################
  
  OResRep[[Simrep]]=list(simuRes1,simuRes2,simuRes3,simuResIni,I,H,Z,R,simuomega1,simuomega2,simuomega3,tildeomegaIni,Isave,Missingrate)
  
}


#Sys.time()-timerep

saveRDS(OResRep,file=paste0('HosCovidCase4', '_irep', irep, '.RDS'))
