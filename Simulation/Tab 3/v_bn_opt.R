v_bn_optim<-function(batch = batch){
  
  
  bn_per<-c()
  md <- as.matrix(dist(batch)^2)
  
  for(i in 1:nrow(batch)){
    
    group_id_aux <- rep(0,nrow(batch))
    group_id_aux[i]<-1
    
    bn_per[i]<-bn(group_id_aux, md = md)
  }
  
  return(var(bn_per))
  
  
}