library(forecast)
library(tsqn)
library(uclust)
source("v_chart_funs.R")
source("fun_heteroc_sim.R")


#### default for "fun_heteroc_sim.R"
# n=1000,I=20,Inew=500,sigt0=1,sigt.new=10,LagMax=200
# confidence.level=0.95,phi=0.4,phi.new=0.4,print.all=FALSE
Rep=2

Iv <- c(20)
nv <- c(500)
Inew <- 300

phi=0.1;phi.new=0.1

sigt0=1;sigt.new=1 #factor of variance 
# t=1:n fator_t.new=(sigt.new/n)*t

perc_T2 <- matrix(0,ncol=length(Iv),nrow=length(nv))



colnames(perc_T2) <- Iv
rownames(perc_T2) <- nv

perc_Va <- perc_T2;  perc_Vd <- perc_T2
perc_Va_sd <- perc_T2;  perc_Vd_sd <- perc_T2
perc_T2_sd <- perc_T2
ARL_T2<-perc_T2

ARL_Va <- perc_T2;  ARL_Vd <- ARL_T2
ARL_Va_sd <- ARL_T2;  ARL_Vd_sd <- ARL_T2
ARL_T2_sd <- ARL_T2

for(i in 1:length(nv)){
  n <- nv[i]
  for (j in 1:length(Iv)) {
    I <- Iv[j]
    p_T2 <- c();  p_Va <- c();  p_Vd <- c()
    
    A_T2<-c(); A_Va<-c(); A_Vd<-c()
    
    for(r in 1:Rep){
      res <- sim_fun_heterok(n=n,I=I,Inew=Inew,phi=phi,phi.new=phi.new,sigt0=sigt0,sigt.new=sigt.new)
      p_T2[r] <- res$percent_T2
      p_Va[r] <- res$percent_Va
      p_Vd[r] <- res$percent_Vd
      A_T2[r]<-ifelse(res$percent_T2==0,NA, 1/res$percent_T2)
      A_Va[r]<-ifelse(res$percent_Va==0,NA, 1/res$percent_Va)
      A_Vd[r]<-ifelse(res$percent_Vd==0,NA, 1/res$percent_Vd)
    }
    
    perc_T2[i,j] <- mean(p_T2)
    perc_T2_sd[i,j] <- sd(p_T2)
    perc_Va[i,j] <- mean(p_Va)
    perc_Va_sd[i,j] <- sd(p_Va)
    perc_Vd[i,j] <- mean(p_Vd)
    perc_Vd_sd[i,j] <- sd(p_Vd)
    
    ARL_T2[i,j] <- mean(A_T2, na.rm=T)
    ARL_T2_sd[i,j] <- sd(A_T2, na.rm=T)
    ARL_Va[i,j] <- mean(A_Va, na.rm=T)
    ARL_Va_sd[i,j] <- sd(A_Va, na.rm=T)
    ARL_Vd[i,j] <- mean(A_Vd, na.rm=T)
    ARL_Vd_sd[i,j] <- sd(A_Vd, na.rm=T)
    
    #save.image(file = "sigt_1_phi_01.RData")
  }
}
