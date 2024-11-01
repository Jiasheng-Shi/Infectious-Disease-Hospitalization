library(mvtnorm)
library(EpiEstim)
library(mc2d)
library(ggplot2)
library(nlme)
library(numDeriv)
library(TTR)
library(parallel)

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
## Overall Parameter Setting 

{
  ## length of the Time series/number of observations
  T=120
  ## Generate Covariants
  ## Z is the covariates, Z[,1:(NoCov=2)] assumed to be time variate
  NoCov=2
  
  ##################################
  ## True parameter value:
  ## Oracle beta
  OracleBeta=as.matrix(c(-0.02,-0.125),nrow=NoCov,ncol=1)
  
  ## Oracle phi
  OraclePhi=c(0.7,0.5)
  OraclePhi_0=OraclePhi[1]
  OraclePhi_1=OraclePhi[2]
  
}


##################################
## Load Sampled Data
{
  
  #setwd()
  SampleData=readRDS(file = "sampled_data_case1.rds")
  I_0=SampleData[[1]]
  R_0=SampleData[[2]]
  Z=SampleData[[3]]
  I=SampleData[[4]]
  R=SampleData[[5]]
  H=SampleData[[6]]
  Isave=SampleData[[7]]
  
}

##################################
## Seeking for initial values
{
  ###########################################################################
  
  ### Composite for omega
  ### Fixed I and H
  
  V3compOmega<-function(paraT,omega,I,Z){
    Rfi<-c()
    Lambdafi<-c()
    NumI=c()
    DenI=c()
    
    ## Rfi(R-function-intermediate) and Lambdafi
    Lambdafi<-WMA(c(rep(0,(length(omega)-1)),I_0,c(I)),
                  n=length(omega),
                  wts = rev(omega))[length(omega):(length(omega)+T-1)]
    
    Rfi[1]=paraT[1] + (paraT[2]*log(R_0)) + (Z[1,]%*%paraT[3:(2+ncol(Z))])
    for (r in 2:T) {
      Rfi[r]= (paraT[1]) + (paraT[2]* Rfi[r-1] )+ (Z[r,]%*%paraT[3:(2+ncol(Z))] )  
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
    
    return(result)
  }
  
  V3compOmega2<-function(x){
    paraT=x
    result=V3compOmega(paraT,omega=omega,I=I,Z=Z)
    return(result)
  }
  
  c(OraclePhi,OracleBeta)
  
  timeS=Sys.time()
  simuResIni=optim(fn=V3compOmega2,par = c(0.6,0.5,-0.01,-0.12), 
                lower = c(-1,-1,-1,-1),upper = c(1,1,1,1),
                method = "L-BFGS-B",
                #control = list(pgtol=0 )
  )
  Sys.time()-timeS
  
  simuResIni$par

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
      successP=exp(successP)
      successP=successP/( 1.5* max(successP) )
      successP[successP>1]=1
      
      successPsi=mapply(FUN = rbinom, n=1,size=1, prob=successP)
      
      if (sum(successPsi)<=(TN0/10)){
        successP=successP*(1/max(successP))
        successP[successP>1]=1
        successPsi=mapply(FUN = rbinom, n=1,size=1, prob=successP)
      }
      
      
      hfi[[r]]=cbind(hfi[[r]],successPsi)
      hfi[[r]]=hfi[[r]][hfi[[r]][,length(tildeomegaInput)]==1,]
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
    cellQLV=sum(cellQoverall/cellQN)
    
    resCompLike=list(cellQLV=cellQLV,cellQoverall=cellQoverall,cellQN=cellQN,cellQlog=cellQlog)
    return(resCompLike)
    
  }
  
}

########## Negative CompLike function with fixed omega and tildeomega
########## First Iteration, arbitrary initial values

{

  NCLFOmega1=function(paraTOutput){
    resNCL=-CompLike(paraTOutput,paraT=c(0.1,0.6,-0.25,-0.25),
                     I=I,
                     H=H,
                     Z=Z,
                     omegaOutput=omega,
                     tildeomegaOutput=tildeomega,
                     omegaInput=omega,
                     tildeomegaInput=tildeomega)$cellQLV
    #resNCL=-sum(resNCL$cellQoverall/resNCL$cellQN)
    return(resNCL)
  }
  
  c(OraclePhi,OracleBeta)
  
  timeS=Sys.time()
  simuRes1=nlm(NCLFOmega1,p= c(0.1,0.6,-0.25,-0.25),steptol = 1e-5)
  Sys.time()-timeS

}

########## Negative CompLike function with fixed omega and tildeomega
########## Second Iteration, arbitrary initial values

{

  NCLFOmega2=function(paraTOutput){
    resNCL=-CompLike(paraTOutput,paraT=simuRes1$estimate,
                     I=I,
                     H=H,
                     Z=Z,
                     omegaOutput=omega,
                     tildeomegaOutput=tildeomega,
                     omegaInput=omega,
                     tildeomegaInput=tildeomega)$cellQLV
    #resNCL=-sum(resNCL$cellQoverall/resNCL$cellQN)
    return(resNCL)
  }
  
  c(OraclePhi,OracleBeta)
  
  timeS=Sys.time()
  simuRes2=nlm(NCLFOmega2,p= simuRes1$estimate,steptol = 1e-6)
  Sys.time()-timeS

}

########## Negative CompLike function with fixed omega and tildeomega
########## Third Iteration, arbitrary initial values

{

  NCLFOmega3=function(paraTOutput){
    resNCL=-CompLike(paraTOutput,paraT=simuRes2$estimate,
                     I=I,
                     H=H,
                     Z=Z,
                     omegaOutput=omega,
                     tildeomegaOutput=tildeomega,
                     omegaInput=omega,
                     tildeomegaInput=tildeomega)$cellQLV
    #resNCL=-sum(resNCL$cellQoverall/resNCL$cellQN)
    return(resNCL)
  }
  
  c(OraclePhi,OracleBeta)
  
  timeS=Sys.time()
  simuRes3=nlm(NCLFOmega3,p= simuRes2$estimate,steptol = 1e-6)
  Sys.time()-timeS

}

##################################
## Calculate the estimated instantaneous reproduction number R

{
  
  Rest1=c()
  paraest1=simuRes1$estimate
  Rest1[1]=paraest1[1] + (paraest1[2]*log(R_0)) + (Z[1,]%*%paraest1[3:(2+ncol(Z))])
  for (r in 2:T) {
    Rest1[r]= (paraest1[1]) + (paraest1[2]* Rest1[r-1] )+ (Z[r,]%*%paraest1[3:(2+ncol(Z))] )  
  }
  Rest1=exp(Rest1)
  
  Rest2=c()
  paraest2=simuRes2$estimate
  Rest2[1]=paraest2[1] + (paraest2[2]*log(R_0)) + (Z[1,]%*%paraest2[3:(2+ncol(Z))])
  for (r in 2:T) {
    Rest2[r]= (paraest2[1]) + (paraest2[2]* Rest2[r-1] )+ (Z[r,]%*%paraest2[3:(2+ncol(Z))] )  
  }
  Rest2=exp(Rest2)
  
}

##################################
## Calculate the estimated instantaneous reproduction number R based on Cori

{
  
  Icori=as.vector(I)
  t_start=c(2:118)
  t_end=t_start+2
  si_default=matrix(data=c(0,omega),nrow = (1+length(omega)),ncol = 1 ,byrow = TRUE)
  RestCori=estimate_R(
    Icori,
    method = "si_from_sample",
    si_sample = si_default,
    config = make_config(list(
       t_start = t_start, 
       t_end = t_end))
  )
  
}

##################################
## ggplot

{
  
  plot(x=c(1:120),y=Rest1,type="l",lty=4)
  lines(x=c(1:120),y=Rest2,lty=2,col=3)
  lines(x=c(1:120),y=R,lty=1,col=2)
  lines(x=c(7:119),y=RestCori$R$`Mean(R)`[5:117],col=4)
  
  sum( (RestCori$R$`Mean(R)`[1:113] - R[7:119])^2  )
  
  
  estRplot=data.frame(diff=c( (Rest1-R)[7:119], RestCori$R$`Mean(R)`[5:117]-R[7:119] ) , 
                      time= c( c(7:119), c(7:119) ),
                      cluster= c( rep("Proposed Method",113),rep("Cori's Method",113) ) )
  
  
  barcharts1<- ggplot(estRplot, aes(x = time , y= diff, fill = cluster)) +
    geom_bar(position="dodge", stat = "identity")
  barcharts1
  
  estRplot=data.frame("Estimated Bias"=c( (Rest1-R)[46:119],(Rest2-R)[46:119], 
                              RestCori$R$`Mean(R)`[44:117]-R[46:119]  ) , 
                      Days= c( c(46:119), c(46:119), c(46:119) ),
                      "Estimation Method"= c( rep("Proposed Method",74), 
                                  rep("Proposed Method second iteration",74),
                                  rep("Cori's Method",74) ) )
  
  barcharts2<-estRplot%>%filter(Estimation.Method != "Proposed Method second iteration")%>%
    ggplot(aes(x = Days , y= Estimated.Bias, fill = Estimation.Method)) +
    geom_bar(position="dodge", stat = "identity")+
    theme(legend.position="top")+
    #scale_fill_manual(values=c("#EDA6C4","#AE98AA","#92A1CF"))+
    #scale_fill_manual(values=c("#EEA6B7","#BCCF90","#51C4D3"))+
    #scale_fill_manual(values=c("#894E64","#BE7E4A","#B14B28"))+
    #scale_fill_manual(values=c("#EDA6C4","#615F74","#9E725A"))+
    #scale_fill_manual(values=c("#EDA6C4","#AE98AA","#656255"))+
    #scale_fill_manual(values=c("red","blue","blue"))+
    labs(x="Days", y="Estimation Bias",fill = "Estimation Method: ")
  #coord_cartesian(ylim=c(-0.025,0.025))
  
  barcharts2
  
}

#saveRDS(list(simuRes1,simuRes2,simuResIni,I,H,Z,barcharts1,barcharts2), file='HosCovidCase1.RDS')

