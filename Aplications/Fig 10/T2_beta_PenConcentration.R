## Data downloaded from 
#https://data.mendeley.com/datasets/pdnjz7zz5x/1   (before 05/2022)

library(latex2exp)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(tseries)
library(FitAR)
#function to plot
mean_cl_quantile <- function(x, q = c(0.1, 0.9), na.rm = TRUE){
  dat <- data.frame(y = mean(x, na.rm = na.rm),
                    ymin = quantile(x, probs = q[1], na.rm = na.rm),
                    ymax = quantile(x, probs = q[2], na.rm = na.rm))
  return(dat)
}
file1<-"100_Batches_IndPenSim_V3.csv"
# peni2 <- read.csv(file1,stringsAsFactors=FALSE, header=T, nrows=5) 
# 
# colnames(peni2)[38]
# dim(peni2)
# 
# btc_ref <- fread(file1, select = c(36)) #1 to 100
# btc_ref_control<-fread(file1, select = c(37))# 0 and 1's
# timeh<-fread(file1, select = 1)
# time<-timeh*10/2
# 
# #Batches 1-30: Controlled by recipe driven approach 
# #Batches 31-60: Controlled by operators 
# #Batches 61:90: Controlled by an Advanced Process Control (APC)  solution using the Raman spectroscopy
# #Batches 91:100: Contain faults resulting in process deviations.
# 
# col_btc=14
# 
# df<-data.frame(c(fread(file1, select = col_btc),btc_ref,btc_ref_control,time))
# df[,5]<-c(0,diff(df[,1]))  
# colnames(df)=c("var","ref","group","time","diff_var")
# dfs<-df
# 
# save.image(file="pen_concentration.RData")  
#### loading filtered data

load(file="pen_concentration.RData")
d_ref=61:90
df<-filter(dfs, ref %in% c(d_ref,91:100))
#table(df$ref)
d_inc<-length(d_ref)

p1<-ggplot(df, 
           aes(x=time, y=var, group=ref,  
               color=factor(group))) +
  geom_point(aes(shape=factor(group)), size=0.1) +
  geom_line()
p1

t_max=835

serie_media=c()
for(t in 1:t_max){
  serie_media[t]=mean(df$var[which(df$time==t)][1:d_inc])
}
plot.ts(serie_media)

df2<-df[which(df$time<=t_max),]
df2$serie_media<-c(rep(serie_media,(d_inc+10)))
df2$diff_var<-df2$var-df2$serie_media # 

#df2<-df2[which(df2$time>2),]
df2<-df2[which(df2$time>2),]
df2$ref2<-df2$ref
df2$ref2[which(df2$group==1)]<-"Faulty group"
df2$ref2[which(df2$group==0)]<-"In-control group"
p2<-ggplot(df2, 
           aes(x=time, y=diff_var, group=ref,  
               color=factor(ref2))) +
  geom_point(aes(shape=factor(group)), size=0.001) +
  geom_line()+
  theme(legend.text = element_text(size=12),legend.position = "bottom")+ labs( y="Oxygen Uptake Rate",colour="")

p2<-p2 + guides(fill = FALSE, linetype = FALSE, shape = FALSE)
p2


#### fit
adt=c()
k=0
for( s in d_ref){
  k=k+1
  s1=df2$diff_var[which(df2$ref==s)]
  adt[k]=adf.test(s1)$p.value
}
sum(adt<0.1)# % of time series that don´t have unit root (adf test)



#########################################
d1=d_inc
id1=d_ref
id0=91:100 #### to monitor

po=1:2
coefs_ar=matrix(0,ncol=(length(po)+1),nrow=(d_inc+10))
f1=list()

k=0
for(s in c(d_ref,91:100)){
  k=k+1
  print(k)
  s1=df2$diff_var[which(df2$ref==s)]
  if(adf.test(s1)$p.value>0.05){
    s1<-diff(s1)
  f1[[k]]=arima(s1,order=c(max(po-1),0,0),include.mean = TRUE)
  coefs_ar[k,]=c(1,f1[[k]]$coef)
  }
  f1[[k]]=arima(s1,order=c(max(po),0,0),include.mean = TRUE)
  coefs_ar[k,]=f1[[k]]$coef
}


############## ACF Residuals

dt=list()
boxtest=c()
shaptest=c()
for(j in 1:d1){
  Acf=acf(f1[[j]]$res,plot=FALSE,lag.max = 60)$acf[-1]
  y=Acf;x=1:60
  boxtest[j]=Box.test(f1[[j]]$res)$p.value
  shaptest[j]=shapiro.test(f1[[j]]$res)$p.value
}

sum(shaptest >= 0.05)/d1
sum(boxtest>=0.05)/d1 #H0: The data are independently distributed 



###################
## individual t-test
t_tets_ind=matrix(0,ncol=(length(po)+1),nrow=d1)
for(j in 1:d1){
  t_tets_ind[j,]=f1[[j]]$coef/sqrt(diag(f1[[j]]$var.coef))
}

ap=apply(t_tets_ind[1:d1,],2,function(x){abs(x)>=qnorm(0.975)*1})
colMeans(ap) # criteria : + more than 70% of times were significant in the individual tests




############## plot all coefs  AR(po)
reff=c()
reff[1:d1]="Faulty group"
reff[31:40]="In-control group"
dt=list()
coefs_ar2=coefs_ar[,c(1,2)]


for(j in 1:(d_inc+10)){
  dt[[j]] <- data.table(Order = seq(1,length(po)+1),coefs = coefs_ar2[j,],ref=rep(reff[j],length(po)+1))
}

dt <- rbindlist(dt,idcol = 'Run') 
p1=ggplot(dt, aes(Order, coefs,group=factor(ref),colour=factor(ref))) +
  stat_summary(geom = "line", fun.y = mean) +
  geom_smooth(stat = 'summary', fun.data = mean_cl_quantile,alpha = 0.2,size=0.25,aes(fill = factor(ref)))
p1=p1 +theme(legend.text = element_text(size=12),legend.position = "bottom")+ labs( y="AR coefficients",fill="",group="",colour="")

Par=p1+scale_x_continuous(n.breaks = 3,labels=c("0","1","2"))#tamanho 1100x 500

Par
grid.arrange(p2, Par, nrow = 1, ncol = 2)#tamanho 15 x 8


## T2 hotteling
coef_use<-c(1,2,3)

mean_coefs=apply(coefs_ar[1:d1,coef_use],2,mean)
cov_coefs=cov(coefs_ar[1:d1,coef_use])


T2=c()

for(i in 1:30){
  T2[i] = mahalanobis(
    coefs_ar[i,coef_use],
    center = mean_coefs,
    cov =cov_coefs,
    inverted = FALSE)
}


T2.new=c()

for(i in 1:10){
  T2.new[i] = mahalanobis(
    coefs_ar[i+d_inc,coef_use],
    center = mean_coefs,
    cov =cov_coefs,
    inverted = FALSE)
  
}

summary(T2.new)

IC=30
p=length(coef_use)
T2_Lim_Teo= ((IC+1)/IC)*(p*(IC-1)/(IC-p))*qf(0.99,p,(IC-p))
sum(T2.new>T2_Lim_Teo)

##############
# individuals
mth=matrix(0,nrow=10,ncol=length(coef_use))

for(i in 1:10){
  mth[i,] = (coefs_ar[i+30,coef_use]-mean_coefs)/sqrt(diag(cov_coefs))
}

t_lim<-sqrt(31/30)*(qt(0.99,29))

colMeans( ((mth>t_lim)|(mth<= -t_lim))*1)



##############################Charts

#######################T2_beta 
T2_g=c(T2,T2.new)



vc=(T2_g>T2_Lim_Teo)*1

IC.new=10


df=data.frame(i=1:(IC+IC.new),T2_g=T2_g, vc=vc)

go<-ggplot(df, aes(x=i,y=T2_g))+geom_point()
#go<-go+ggtitle("T2_Beta scores of Dynamic ARMA coefficients")
#go<-go+geom_hline(yintercept = 0)
go<-go+geom_vline(xintercept = IC,linetype="dashed",size=0.6)
go<-go+xlab("batches")+ylab(TeX("$T^2_{\\beta}$ scores"))
go<-go+geom_hline(yintercept = T2_Lim_Teo)
go<-go+geom_point(aes(color=as.factor(vc)),y=df[,2],x=df[,1])+ylim(c(0,max(max(df$T2_g),T2_Lim_Teo)+5))
go<-go+guides(color = FALSE)
go<-go+scale_color_manual(values=c("blue","red"))
g=go  
# export pdf landscape 4 x 9
g
