


sim_fun_heterok=function(n=1000,I=20,Inew=500,sigt0=1,sigt.new=10,LagMax=200,confidence.level=0.95,phi=0.4,phi.new=0.4,print.all=FALSE){
  print("####################################")
  print(paste("n=",n," - - I=",I," - - Inew=",Inew,sep=""))
  
  
  if(print.all==TRUE){
    print("sigt0=");print(sigt0)
    print("sigt.new=");print(sigt.new)
    print("LagMax=");print(LagMax)
    print("confidence.level=");print(confidence.level)
    print("phi=");print(phi)
    print("phi.new=");print(phi.new)
  }
  
  t=1:n
  if(sigt0==1){fator_t=1}else{fator_t=(sigt0/n)*t}
  if(sigt.new==1){fator_t.new=1}else{fator_t.new=(sigt.new/n)*t}
  
  
  X=matrix(0,ncol=n,nrow=I)
  
  l.phi=c()
  l.theta=c()
  
  
  inv.coef.cov="erro"
  while(inv.coef.cov[1]=="erro"){
    
    for(i in 1:I){
      #x <-  c(arima.sim(model=list(ma=theta),n=n/2,sd=1),
      #      arima.sim(model=list(ar=phi),n=n/2,sd=1))
      innov=rnorm(n)*fator_t
      x <-  arima.sim(model=list(ar=phi),n=n,sd=1,innov=innov)
      X[i,] <-x
      arm=auto.arima(x)
      l.phi[i]=length(arm$model$phi)
      #print(arm$model$phi)
      l.theta[i]=length(arm$model$theta)
    }
    order.p=max(l.phi)
    order.q=max(l.theta)
    fit.ar=matrix(0,ncol=(order.p+order.q),nrow=I)
    fit.ar.var=list()
    
    for(i in 1:I){
      #print(i)
      fita="erro"
      while(fita[1]=="erro"){
        tryCatch(
          {fita=arima(X[i,], order=c(order.p,0,order.q),include.mean = FALSE)},
          error = function(e){fita="erro"}
        )
        if(fita=="erro"){
          innov=rnorm(n)*fator_t
          X[i,]<- arima.sim(model=list(ar=phi),n=n,sd=1,innov=innov)
        }
      }
      
      
      
      
      fit.ar[i,] <- fita$coef
      fit.ar.var[[i]] <- fita$var.coef
    }
    coef.mean = apply(fit.ar,2,mean)
    #coef.cov = Reduce('+',fit.ar.var)/IC
    coef.cov = cov(fit.ar)
    tryCatch(
      {inv.coef.cov=solve(coef.cov)},
      error = function(e){inv.coef.cov="erro"}
    )
  }
  
  
  
  
  
  
  ############# Fase II ############
  
  fit.ar.new=matrix(0,ncol=(order.p+order.q),nrow=Inew)
  Xnew=matrix(0,ncol=n,nrow=Inew)
  for(i in 1:Inew){
    
    fita="erro"
    while(fita[1]=="erro"){
      innov.new=rnorm(n)*fator_t.new
      xnew <- arima.sim(model=list(ar=phi.new),n=n,sd=1,innov=innov.new)
      
      tryCatch(
        {fita=arima(xnew, order=c(order.p,0,order.q),include.mean = FALSE)},
        error = function(e){fita="erro"}
      )
    }
    #print(i)
    Xnew[i,]<-xnew
    
    fit.ar.new[i,] <- fita$coef
  }
  
  
  T2.new=c()
  for(i in 1:Inew){
    T2.new[i] = mahalanobis(
      fit.ar.new[i, ],
      center = coef.mean,
      cov =coef.cov,
      inverted = FALSE)
  }
  order.p+order.q->vv
  T2_Lim= ((I+1)/I)*(vv*(I-1)/(I-vv))*qf(confidence.level,vv,(I-vv))
  
  percent_T2=mean(T2.new>T2_Lim)
  
  ##########V chart
  vc=vchart(LagMax=LagMax,X=X,Xnew=Xnew,conf.level=confidence.level)
  percent_Va=vc$r_Va
  percent_Vd=vc$r_Vd
  
  return(list(n=n,I=I,Inew=Inew,percent_T2=percent_T2,percent_Va=percent_Va,percent_Vd=percent_Vd))
}