# load packages
library(pracma)
library(ggplot2)
library(EpiEstim)
library(nlme)
library(dlnm)
library(tsModel)
library(corrplot)
library(mvtnorm)
library(stats)
library(gridExtra)
library(ggridges)
library(lubridate)

# Overall Parameters
R_0=2.5
I_0=200
## transmissibility rate of the epidemic (infectiousness profile)

## tau_0, pre-specified time point where before date tau_0, the MLE will be applied.
tau_0=5 #default
## bias_corr_const, the bias correction constant, default is bias_corr_const=1
bias_corr_const=exp(-0.001/2)

# Functions
## logistic function
logit<-function(x){
  result=log(x/(1-x))
  return(result)
}
##
greater.or.na <- function(myobj, myconst, tocat){
  if(is.na(myobj)) cat(myobj, tocat, '\n')
  return(is.na(myobj) | myobj > myconst)
}


############ Main Function ##########################

QSOEID<-function(Z,I){

      ##################################################

      ### Lambda[t]
      Lambda<-matrix(data = NA, nrow=1,ncol=T)
      Lambda[1]=I_0*Omega[1]
      for (t in 2:25) { Lambda[t]=(I_0*Omega[t])+ (I[1:(t-1)]%*%Omega[(t-1):1] ) }
      for (t in 26:T) { Lambda[t]= I[(t-1):(t-25)]%*%Omega[1:25]   }

      ##################################################

      # Estimating EstR[1:tau_0] based on MLE
      EstR<-matrix(data = NA, nrow=1,ncol=T)
      EstR[1]=max(1,(I[1]/(I_0*Omega[1]) ) )
      for (t in 2:tau_0) { EstR[t]=max(1,  (I[t]/( (I_0*Omega[t])+ (I[1:(t-1)]%*%Omega[(t-1):1] ) ) ) )  }

      #################################################
      
      ### barZ[t], arverage of Z[i,]'s from date 1 to date t.
      barZ<-Z
      for (t in 2:T) { barZ[t,]=colSums(Z[1:t,])/t }
      
      ### WTilde, weight  \Big( \sum_{i=1}^t (Z_i-\bar{Z})(Z_i-\bar{Z})^T \Big)^{-1}
      WTilde<-list()
      WTilde[[1]]<-matrix(data=0, nrow=NoCov, ncol=NoCov)
      for (t in 2:tau_0) {
        WTilde[[t]]<-matrix(data=0, nrow=NoCov, ncol=NoCov)
        for (i in 1:t) {
          WTilde[[t]]=WTilde[[t]]+(Z[i,]-barZ[t,])%*%t((Z[i,]-barZ[t,]))
        }
        WTilde[[t]]=solve( WTilde[[t]] +diag(c(rep(0.1,NoCov)),nrow=NoCov ))
      }
      for (t in 3:T) {
        WTilde[[t]]<-matrix(data=0, nrow=NoCov, ncol=NoCov)
        for (i in 1:t) {
          WTilde[[t]]=WTilde[[t]]+(Z[i,]-barZ[t,])%*%t((Z[i,]-barZ[t,]))
        }
        WTilde[[t]]=solve( WTilde[[t]] )
      }

      ###
      YTilde<-matrix(data = NA, nrow=1,ncol=T)
      for (t in 1:tau_0) { YTilde[t]=log(EstR[t])}
      BetaTilde<-matrix(data=NA, nrow=T, ncol=NoCov)
      ZYTilde<-matrix(data = NA, nrow=T,ncol=NoCov)
      ### EstPhi, EstBeta, for instore the estimated regression parameter
      EstPhi<-matrix(data = NA,nrow=T,ncol = 2)
      EstBeta<-matrix(data = NA, nrow=T,ncol = NoCov)
      ### Intermediate variables, ZYHat
      ZYHat=ZYTilde
      
      ### Define the profile likelihood
      ell<- function(phi,k){
        ZYTilde[(k-1),]=( log(EstR[,1]) -phi[2]*log(R_0) -phi[1]  )%*%(Z[1,]-barZ[(k-1),])
        EstR[,1]=exp(phi[1]+phi[2]*log(R_0))
        for (i in 2:(k-1)) {
          ZYTilde[(k-1),]=ZYTilde[(k-1),]+ (log(EstR[,i])-phi[2]*log(EstR[,(i-1)]) - phi[1]  ) %*%(Z[i,]-barZ[(k-1),])
        }
        BetaTilde[k,]=ZYTilde[(k-1),]%*%WTilde[[(k-1)]]
        for (i in (tau_0+1):k){
          YTilde[,i]=phi[1]+phi[2]*log(EstR[,(i-1)])+t(Z[i,])%*%BetaTilde[k,]
        }
        result=0
        for (j in (tau_0+1):k) {
          result=result+I[j]*YTilde[,j]-bias_corr_const*(exp(YTilde[,j]))*Lambda[j]
        }
        return(-result)
      }
 
      ### Minimize over the minus profile log-likelihood
      for (t in (tau_0+1):T) {
        EstPhi[t,]=nlminb(c(0.05,0.7),ell,k=t,lower = c(-5,0.3),upper = c(5,0.95))$par
        
        ZYHat[(t-1),]=( log(EstR[,1]) -EstPhi[t,2]*log(R_0) -EstPhi[t,1]  )%*%(Z[1,]-barZ[(t-1),])
        for (i in 2:(t-1)) {
          ZYHat[(t-1),]=ZYHat[(t-1),]+ (log(EstR[,i])-EstPhi[t,2]*log(EstR[,(i-1)]) - EstPhi[t,1]  ) %*%(Z[i,]-barZ[(t-1),])
        }
        EstBeta[t,]=ZYHat[(t-1),]%*%WTilde[[(t-1)]]
        
        EstR[,t]=exp( EstPhi[t,1]+EstPhi[t,2]*log(EstR[,(t-1)])+t(Z[t,])%*%EstBeta[t,] ) 
        
        if(greater.or.na(abs(EstR[,t]-EstR[,(t-1)]), 5, 'L182?')){
          EstPhi[t,]=EstPhi[(t-1),]
          EstBeta[t,]=EstBeta[(t-1),]
          EstR[,t]=EstR[,(t-1)]
        }
      }	
      
    
      ### Calculate the bootstrap confidence interval for parameters.
      
      restPara<-matrix(data = NA,nrow=1,ncol = 4)
      restPara[1,1]=EstPhi[T,1]
      restPara[1,2]=EstPhi[T,2]
      restPara[1,3]=EstBeta[T,1]
      restPara[1,4]=EstBeta[T,2]
      
      resQSOEID=list(restPara=restPara,EstR=EstR)
      return(resQSOEID)
}












