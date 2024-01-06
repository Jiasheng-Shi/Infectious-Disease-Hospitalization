# load packages
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

#####Model 1, corrctly specified model##########################################################

###################################
## Generate the instantaneous reproduction number 

{
  set.seed(10)
  ## length of the Time series/number of observations
  T=120
  ## Generate Covariants
  ## Z is the covariates, Z[,1:(NoCov=2)] assumed to be time variate
  NoCov=2
  
  Z<-matrix(data=NA, nrow=T, ncol=NoCov)
  
  for (t in 1:T) {Z[t,1]=5-(T/8)+((2*t)/8)+rnorm(1,mean=0,sd=2.5) } 
  
  logit=function(p){return(log(p/(1-p)))}
  Z[,2]=logit( runif(T,min=0.01,max=0.99) )+2 
  
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
  set.seed(10)
  I_0=50
  
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
  set.seed(10)
  h=rmultinomial(n=T,size = I, prob = tildeomega)
  H=matrix(data=NA,nrow=1,ncol=T)
  for (i in 1:T) {
    H[i]=0
    for (j in max(1,i+2-length(tildeomega)):i ) {
      H[i]=H[i]+h[j,(i+1-j)]
    }
  }
  
}

SampleData=list(I_0,R_0,Z,I,R,H,Isave)

saveRDS(SampleData,file="sampled_data_case1.rds")













