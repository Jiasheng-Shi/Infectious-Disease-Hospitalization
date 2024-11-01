library(mvtnorm)
library(EpiEstim)
library(mc2d)
#library(ggplot2)
library(nlme)
#library(numDeriv)
library(TTR)
#library(parallel)
library(dplyr)

###############################
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


####################################
### For Table
{
  
  OResRep=list()
  OResRep[[1]]=matrix(data=NA, nrow = 1000,ncol = 4)
  OResRep[[2]]=matrix(data=NA, nrow = 1000,ncol = 24)
  OResRep[[3]]=matrix(data=NA, nrow = 1000,ncol = 7)
  OResRep[[4]]=matrix(data=NA, nrow = 1000,ncol = 120)
  OResRep[[5]]=matrix(data=NA, nrow = 1000,ncol = 120)
  OResRep[[6]]=matrix(data=NA, nrow = 1000,ncol = 120)
  OResRep[[7]]=matrix(data=NA, nrow = 1000,ncol = 120)
  OResRep[[8]]=list()
  
  for (Simrep in 1:10) {
    interaa=readRDS(file = paste0('HosCovidCase4', '_irep', Simrep, '.RDS'))
    for (irep in 1:100){
      OResRep[[1]][(Simrep*100+irep-100),]=interaa[[irep]][[1]]$estimate #para
      OResRep[[2]][(Simrep*100+irep-100),]=interaa[[irep]][[9]]$par #omega
      
      OResRep[[3]][(Simrep*100+irep-100),]=interaa[[irep]][[12]]   #tildeomgea
      OResRep[[4]][(Simrep*100+irep-100),]=interaa[[irep]][[13]]   #Isave
      OResRep[[5]][(Simrep*100+irep-100),]=interaa[[irep]][[14]]   #missingrate
      OResRep[[6]][(Simrep*100+irep-100),]=interaa[[irep]][[8]] #R
      OResRep[[7]][(Simrep*100+irep-100),]=interaa[[irep]][[5]]   #I with missing
      OResRep[[8]][[(Simrep*100+irep-100)]]=matrix(data=NA, nrow = 120,ncol = 2)
      OResRep[[8]][[(Simrep*100+irep-100)]]=interaa[[irep]][[7]]   #Z
      
    }
    
  }
  
  OResRep[[2]]=cbind(OResRep[[2]],1-rowSums(OResRep[[2]]))
  
}

###################################
##  the instantaneous reproduction number 
## HosPaper
## R[t], the instantaneous reproduction number at time t.
{ 
   
  R1O<-matrix(data = NA, nrow=1,ncol=T)
  ### R[0], the instantaneous reproduction number at time 0.
  R_0=2.5
  
  Z1=OResRep[[8]][[1]]
  OraclePhi_0=0.7
  OraclePhi_1=0.5
  OracleBeta=c(-0.02,-0.125)
  T=120
  #epsilon<-rep(0,T)
  R1O[1]=exp( (OraclePhi_0) + (OraclePhi_1*log(R_0)) + (Z1[1,]%*%OracleBeta)  )
  for (t in 2:T) {
    R1O[t]=exp( (OraclePhi_0) + (OraclePhi_1*log(R1O[(t-1)]))+ (Z1[t,]%*%OracleBeta)  )
  }
  plot(x=c(1:120),y=R1O[1,],type = "l",lty=4)
  lines(x=c(0:120),y=rep(1,121),type="l",lty=4)
  
  R11<-matrix(data = NA, nrow=1,ncol=T)
  #pa1=OResRep[[1]][1,]
  pa1=colMeans(OResRep[[1]])
  R11[1]=exp( (pa1[1]) + (pa1[2]*log(R_0)) + (Z1[1,]%*%pa1[3:4])  )
  for (t in 2:T) {
    R11[t]=exp( (pa1[1]) + (pa1[2]*log(R11[(t-1)]))+ (Z1[t,]%*%pa1[3:4])  )
  }
  lines(x=c(1:120),y=R11,col=2)
  
}


####################################
## Baseline Paper
{
  
  source("QSOEID.R")
  
  NoCov=2
  tunningl=T
  rep=1
  #Z=OResRep[[8]][[1]]
  restrial1=list()
  restrial1[[1]]=matrix(data = NA, nrow = 1000,ncol = 4)
  restrial1[[2]]=matrix(data = NA, nrow = 1000,ncol = 120)
  
  time1=Sys.time()
  for (i in 1:1000) {
    Omega=OResRep[[2]][i,]

    interaaa=QSOEID(Z=OResRep[[8]][[1]],I=OResRep[[7]][i,])
    restrial1[[1]][i,]=interaaa$restPara
    restrial1[[2]][i,]=interaaa$EstR
  }
  Sys.time()-time1
  
  lines(x=c(1:120),y=restrial1[[2]][1,],col=4)
  lines(x=c(1:120),y=colMeans(restrial1[[2]]),col=5)
  
}

#############################################################

BootQSOEID=restrial1

saveRDS(BootQSOEID,file = "BootQSOEID.RDS")


##############################################################
### For EstR pic
{

  restrial1=readRDS(file = "BootQSOEID.RDS")
  
  EstRCase4<-matrix(data=NA, nrow=1000, ncol = T)
  for (i in 1:1000) {
    
    interpC4=OResRep[[1]][i,]
    EstRCase4[i,1]=(interpC4[1]) + (interpC4[2]*log(R_0)) + (Z1[1,]%*%interpC4[3:4])
    for (t in 2:T) {
      EstRCase4[i,t]=(interpC4[1]) + (interpC4[2]*EstRCase4[i,(t-1)]) + (Z1[t,]%*%interpC4[3:4])
    }
    EstRCase4[i,]=exp(EstRCase4[i,])
    
  }
  
  Quantile5F=function(x){
    result=sort(x)[51]
    return(result)
  }
  
  Quantile95F=function(x){
    result=sort(x)[950]
    return(result)
  }
  
  reg.EstR=data.frame(Time=c(1:T))
  #reg.EstR$Oracle=R1O[1,]
  reg.EstR$Oracle=R1O[1,]+rnorm(T, mean = 0,sd=0.05)
  reg.EstR$EstR=colMeans(EstRCase4)
  reg.EstR$QSOEID=colMeans(restrial1[[2]])
  
  reg.EstR$Low=apply(EstRCase4, 2, Quantile5F)
  reg.EstR$Up=apply(EstRCase4, 2, Quantile95F)
  
  reg.EstR$exLow=apply(EstRCase4, 2, min)
  reg.EstR$exUp=apply(EstRCase4, 2, max)
  
  reg.EstR=rbind(reg.EstR,reg.EstR,reg.EstR)
  reg.EstR$group=rep(c("Estimated Rt, incorporating hospital admission data",
                       "Oracle instantaneous reproduction number",
                       "Estimated Rt, without incorporating hospital admission data"),each=120)
  reg.EstR$Oracle[1:120]=reg.EstR$EstR[1:120]
  reg.EstR$Oracle[241:360]=reg.EstR$QSOEID[1:120]
  reg.EstR$EstR[121:240]=reg.EstR$Oracle[121:240]
  reg.EstR$EstR[241:360]=reg.EstR$QSOEID[1:120]
  reg.EstR$QSOEID[1:120]=reg.EstR$EstR[1:120]
  reg.EstR$QSOEID[121:240]=reg.EstR$Oracle[121:240]
  
  library(latex2exp)
  
  ggplot(data=reg.EstR,aes(x=Time))+
    #geom_line(aes(y=Oracle,colour=group,lty=group))+
    #geom_line(aes(y=EstR),col=2)+
    geom_line(aes(y=EstR,colour=group,lty=group))+
    geom_ribbon(aes(ymin=Low,ymax=Up),col="gray",alpha=0.2)+
    xlab("Time (day)")+
    ylab("instantaneous reproduction number")+
    theme(
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.position = c(0.75,0.8)
    )+
    scale_x_continuous(name='Time', breaks=c(30, 60, 90, 120), 
                       labels=c(30, 60, 90, 120), limits=c(15,120) )+
    ylim(0,7)+
    scale_color_manual(values=c("#CC0033", "dodgerblue","black"))+
    scale_linetype_manual(values=c("solid","twodash","dashed"))+
    geom_line(aes(y=1),color="#006699",lty=2)+
    annotate(geom='text', x=20, y=1.2, label=TeX("$R_t = 1$", output='character'), parse=TRUE)
  
  #1045*410
  
}

#######################################
### For Omega pic
{
  
  reg.Omega=data.frame(TSI=c(1:25))
  reg.Omega$Oracle=omega
  reg.Omega$Est=colMeans(OResRep[[2]])
  
  Quantile5F=function(x){
    result=sort(x)[51]
    return(result)
  }
  
  Quantile95F=function(x){
    result=sort(x)[950]
    return(result)
  }
  
  reg.Omega$Low=apply(OResRep[[2]], 2, Quantile5F)
  reg.Omega$Up=apply(OResRep[[2]], 2, Quantile95F)
  
  ggplot(data=reg.Omega,aes(x=TSI))+
    geom_line(aes(y=Oracle))+
    geom_line(aes(y=Est),col=2)+
    geom_ribbon(aes(ymin=Low,ymax=Up),col="gray",alpha=0.2)
  
  reg.Omega$diff=reg.Omega$Est-reg.Omega$Oracle
  reg.Omega$group=sign(reg.Omega$diff)
  
  ggplot(reg.Omega, aes(
    x = TSI, 
    y = diff,
    fill = group) )+
    geom_col() + 
    # 
    guides(fill = "none")
  
  reg.Omega2=OResRep[[2]]-matrix(data = rep(omega,1000),nrow = 1000,ncol = 25,byrow = TRUE)
  reg.Omega2=data.frame(diff=c(reg.Omega2),TSI=rep(c(1:25),each=1000))
  
  ggplot(reg.Omega2, aes(
    group = TSI, 
    y = diff))+
    #fill = hi_lo)) +
    geom_boxplot() + 
    # 
    guides(fill = "none")
  
}



