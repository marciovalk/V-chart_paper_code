ma <- function(arr, n=15){
  res = arr
  
  for(i in (floor(n/2)+1):(length(arr)-floor(n/2))){
    
    res[i] = mean(arr[(i-floor(n/2)):(i+floor(n/2))])
    
  }
  res
}




vchart=function(LagMax=LagMax,X=X,Xnew=Xnew,conf.level=conf.level){
  n=dim(X)[2]
  n1=dim(X)[1]
  x=matrix(0,ncol=n,nrow=n1)
  y=x
  Re=dim(Xnew)[1]
  xx=matrix(0,ncol=n,nrow=Re);yy=xx
  acf_0=matrix(0,nrow=n1,ncol=LagMax)
  acf_1=matrix(0,nrow=Re,ncol=LagMax)
  
  
  for(r in 1:n1){
    
    y[r,]<-X[r,]-ma(X[r,])
    x[r,]<-X[r,]/sd(X[r,])
    acf_0[r,]=acf((y[r,1:n]),lag.max=LagMax,plot=FALSE)$acf[-1]
  }
  x_h0<-sweep(x, 2, colMeans(x))#retira a serie media dos dados
  
  
  for(r in 1:Re){
   
    yy[r,]<-Xnew[r,]-ma(Xnew[r,])
    xx[r,]<-Xnew[r,]/sd(Xnew[r,])
    acf_1[r,]=acf(yy[r,1:n],lag.max=LagMax,plot=FALSE)$acf[-1]
  }
  x_h1<-sweep(xx, 2, colMeans(x))#retira a serie media(h0) dos dados (h1)
  
  d_a<-as.matrix(((dist(rbind(acf_0,acf_1)))^2))
  d_d<-as.matrix(dist(rbind(x_h0,x_h1))^2)
  


###########################
U0=c()
U0d=c()
for(i in 1:n1){
  group_id=rep(1,n1)
  group_id[i]=0
  U0[i]=bn(group_id, md=d_a[1:n1,1:n1])#dinamic
  U0d[i]=bn(group_id, md=d_d[1:n1,1:n1])#drift
}
U1=c()
U1d=c()
for(i in 1:Re){
  group_id=c(rep(1,(n1-1)),0)
  U1[i]=bn(group_id, md=d_a[c(1:(n1-1),(i+n1)),c(1:(n1-1),(i+n1))])
  U1d[i]=bn(group_id, md=d_d[c(1:(n1-1),(i+n1)),c(1:(n1-1),(i+n1))])
}


#########################% fora e ARL

#z=qnorm(conf.level+(1-conf.level)/2)

z=abs(qnorm(conf.level))

s_Va= sum((U1/sd(U0))>z)
va=U1/sd(U0)

s_Vd= sum((U1d/sd(U0d))>z)
vd=U1d/sd(U0d)

r_Va=s_Va/Re

r_Vd=s_Vd/Re

arl_Va=ifelse(s_Va==0,NA,1/(sum(s_Va)/Re)) 

arl_Vd=ifelse(s_Vd==0,NA,1/(sum(s_Vd)/Re)) 

 
return(list(r_Va=r_Va,r_Vd=r_Vd,arl_Va=arl_Va,arl_Vd=arl_Vd,va=va,vd=vd))

}# end
