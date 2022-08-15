
ma <- function(arr, n=10){
  res = arr
  
  for(i in (floor(n/2)+1):(length(arr)-floor(n/2))){
    
    res[i] = mean(arr[(i-floor(n/2)):(i+floor(n/2))])
    
  }
  res
}
library('uclust')
source("v_bn_opt.R")
source("mpca_arma_fun.R")


# num instants
n_col <-2000 #  20 or 2000
# num bat fase I
n_orig <- c(500) #or c(50,100,1000)
# num bat fase 2
n_fase2 <-500
# acf
lag_max_acf<- n_col-1
alfa_h1<-c(0.5)
theta<-alfa_h1  # MA coeficient

REEP=2  #replications times

# disturbance -0.1
lista_das_listas_acf<-list()
lista_das_listas_drift<-list()
lista_das_listas_mpa_T<-list()
lista_das_listas_mpa_Q<-list()

inicio_0<-Sys.time()
for( n_ite_norig in 1:length(n_orig)){
  
  lista_alfa_acf <- list()
  lista_alfa_drift <- list()
  mpa_list_T2<-list()
  mpa_list_Q<-list()
  print("norig")
  print(n_orig[n_ite_norig])
  
  for( alfas_utili in 1:length(alfa_h1)){
    print("alfa")
    print(alfa_h1[alfas_utili])
    
    resultados_replicacao_acf <- c()
    resultados_replicacao_drift <- c()
    arl_Q<-c()
    arl_T2<-c()
    inicio<-Sys.time()
    for(repl in 1:REEP){
      inicio_rep <- Sys.time()
      # base
      # fase 1
      
      mdserh0<-matrix(rep(NA), 
                      nrow = n_orig[n_ite_norig],
                      ncol = n_col)
      
      for(i in 1:n_orig[n_ite_norig] ){
        
        mdserh0[i,] <- arima.sim(n=n_col, model = list(ar= 0.2))
        
      }
      
      
      mdserh1<-matrix(rep(NA), 
                      nrow = n_fase2,
                      ncol = n_col)
      
      for(i in 1:n_fase2){
        
        mdserh1[i,] <- arima.sim(n=n_col,
                                 model = list(ma = theta))
        
      }
      
      #------------------
      # ACF -------------
      #------------------
      
      
      mdh0_acf <- matrix(rep(NA), 
                         nrow = n_orig[n_ite_norig],
                         ncol = lag_max_acf+1)
      
      for(i in 1:n_orig[n_ite_norig]){
        
        serie<-mdserh0[i,]
        serie_ma<-serie-ma(serie)
        aux<-NULL
        aux<-acf(serie_ma, plot = FALSE, lag.max =lag_max_acf)$acf
        
        mdh0_acf[i,]<-aux 
        
      }
      
      mdh1_acf <- matrix(rep(NA), 
                         nrow = n_fase2,
                         ncol = lag_max_acf+1)
      
      for(i in 1:n_fase2){
        
        serie<-mdserh1[i,]
        serie_ma<-serie-ma(serie)
        aux<-NULL
        aux<-acf(serie_ma, plot = FALSE, lag.max =lag_max_acf)$acf
        mdh1_acf[i,]<-aux
      }
      
      #------------------
      # Drift------------
      #------------------
      
      mean_serieh0 <-colMeans(mdserh0)
      mdserh0_drift<-mdserh0[]-mean_serieh0[col(mdserh0)]
      mdserh1_drift<-mdserh1[]-mean_serieh0[col(mdserh1)]
      
      # calc bn --------------------
      
      bn_new_pad_acf<-c()
      bn_new_pad_drift<-c()
      v_bn_calc_acf<-v_bn_optim(mdh0_acf)
      v_bn_calc_drift<-v_bn_optim(mdserh0_drift)
      
      md_full_acf <- as.matrix(dist(rbind(mdh0_acf,mdh1_acf))^2)
      nra<-nrow(mdh0_acf)
      inicio_acf<-Sys.time()
      
      for(i in 1:nrow(mdh1_acf)){
        # bn padronizado
        mdaux<- md_full_acf[c(1:nra,(nra+i)),c(1:nra,(nra+i))]
        
        bn_new_pad_acf[i] <- uclust::bn(group_id = c(rep(0,nra),1),md = mdaux)/sqrt(v_bn_calc_acf)
      }
      
      fim_acf<-Sys.time()
      # 
      
      if(repl < 5){
        print("Acf:")
        print(fim_acf-inicio_acf)
      }
      
      inicio_drift<-Sys.time()
  
      md_full_drift <- as.matrix(dist(rbind(mdserh0_drift,mdserh1_drift))^2)
      nra<-nrow(mdserh0_drift)
      
      
      for(i in 1:nrow(mdserh1_drift)){
        # bn pad
        mdaux<- md_full_drift[c(1:nra,(nra+i)),c(1:nra,(nra+i))]
        bn_new_pad_drift[i] <- uclust::bn(group_id = c(rep(0,nra),1),
                                          md = mdaux)/sqrt(v_bn_calc_drift)
      }
      
#
      fim_drift<-Sys.time()
      if(repl < 5){
        print("Drift:")
        print(fim_drift-inicio_drift)
      }
      resultados_replicacao_acf[repl]<-ifelse(sum(bn_new_pad_acf > -qnorm( p = 0.05) ) == 0,NA, 1/(sum(bn_new_pad_acf > -qnorm( p = 0.05)  )/nrow(mdh1_acf))) 
      resultados_replicacao_drift[repl]<-ifelse(sum(bn_new_pad_drift > -qnorm( p = 0.05) ) ==0,NA, 1/(sum(bn_new_pad_drift > -qnorm( p = 0.05)  )/nrow(mdserh1_drift))) 

      
      inicio_arl<-Sys.time()
      arl_T2[repl]<-mpca_arma_fun(x = mdserh0,
                                  xnew = mdserh1,
                                  I = n_orig[n_ite_norig],
                                  Inew = nrow(mdserh1),
                                  n=ncol(mdserh1),
                                  confidence_alfa = 0.95)$arl_T2
      
      
      arl_Q[repl]<-mpca_arma_fun(x = mdserh0,
                                 xnew = mdserh1,
                                 I = n_orig[n_ite_norig],
                                 Inew = nrow(mdserh1),
                                 n=ncol(mdserh1),
                                 confidence_alfa = 0.95)$arl_Q
      fim_arl<-Sys.time()
      if(repl <5){
        print("Arl")
        print(fim_arl - inicio_arl)
      }
      
      fim_rep <- Sys.time()
      print("time of replication:")
      print(fim_rep-inicio_rep)
    }
    
    print(round(Sys.time()-inicio))
    lista_alfa_acf[[alfas_utili]]<-resultados_replicacao_acf
    lista_alfa_drift[[alfas_utili]]<-resultados_replicacao_drift
    mpa_list_T2[[alfas_utili]]<-arl_T2
    mpa_list_Q[[alfas_utili]]<-arl_Q
    
  }
  lista_das_listas_mpa_T[[n_ite_norig]]<-mpa_list_T2
  lista_das_listas_mpa_Q[[n_ite_norig]]<-mpa_list_Q
  lista_das_listas_acf[[n_ite_norig]]<-lista_alfa_acf
  lista_das_listas_drift[[n_ite_norig]]<-lista_alfa_drift
}
fim<-Sys.time()
fim-inicio_0
# 
# acf ------------------------------------------

for(j in 1:length(lista_das_listas_acf)){
  for( i in 1:length(lista_alfa_acf)){

    if(i==1 & j== 1){
      media <- mean(lista_das_listas_acf[[j]][[i]],na.rm =TRUE)
      desvpad <- sd(lista_das_listas_acf[[j]][[i]],na.rm =TRUE)
      median <- median(lista_das_listas_acf[[j]][[i]],na.rm =TRUE)
      n_origem <- n_orig[[j]]
      coeficiente <- alfa_h1[[i]]
      
      df_resultados_acf<-data.frame(media,
                                    desvpad,
                                    median,
                                    n_origem,
                                    coeficiente)
      
    }else{
      
      aux<-NULL
      media <- mean(lista_das_listas_acf[[j]][[i]],na.rm =TRUE)
      desvpad <- sd(lista_das_listas_acf[[j]][[i]],na.rm =TRUE)
      median <- median(lista_das_listas_acf[[j]][[i]],na.rm =TRUE)
      n_origem <- n_orig[[j]]
      coeficiente <- alfa_h1[[i]]
      
      aux<-data.frame(media,
                      desvpad,
                      median,
                      n_origem,
                      coeficiente)
      df_resultados_acf<-rbind(df_resultados_acf,aux)
    }
 
  }
}




#  drift ------------------------------------------

for(j in 1:length(lista_das_listas_drift)){
  for( i in 1:length(lista_alfa_drift)){
    
    if(i==1 & j== 1){
  
      media <- mean(lista_das_listas_drift[[j]][[i]],na.rm =TRUE)
      desvpad <- sd(lista_das_listas_drift[[j]][[i]],na.rm =TRUE)
      median <- median(lista_das_listas_drift[[j]][[i]],na.rm =TRUE)
      n_origem <- n_orig[[j]]
      coeficiente <- alfa_h1[[i]]
      
      df_resultados_drift<-data.frame(media,
                                      desvpad,
                                      median,
                                      n_origem,
                                      coeficiente)
      
    }else{
      
      aux<-NULL
      media <- mean(lista_das_listas_drift[[j]][[i]],na.rm =TRUE)
      desvpad <- sd(lista_das_listas_drift[[j]][[i]],na.rm =TRUE)
      median <- median(lista_das_listas_drift[[j]][[i]],na.rm =TRUE)
      n_origem <- n_orig[[j]]
      coeficiente <- alfa_h1[[i]]
      
      aux<-data.frame(media,
                      desvpad,
                      median,
                      n_origem,
                      coeficiente)
      df_resultados_drift<-rbind(df_resultados_drift,aux)
      
    }
  }
}



#  alt T ------------------------------------------

for(j in 1:length(lista_das_listas_mpa_T)){
  for( i in 1:length(lista_alfa_drift)){
    if(i==1 & j== 1){
      media <- mean(lista_das_listas_mpa_T[[j]][[i]],na.rm =TRUE)
      desvpad <- sd(lista_das_listas_mpa_T[[j]][[i]],na.rm =TRUE)
      median <- median(lista_das_listas_mpa_T[[j]][[i]],na.rm =TRUE)
      n_origem <- n_orig[[j]]
      coeficiente <- alfa_h1[[i]]
      
      df_resultados_alT<-data.frame(media,
                                    desvpad,
                                    median,
                                    n_origem,
                                    coeficiente)
      
    }else{
      
      aux<-NULL
      media <- mean(lista_das_listas_mpa_T[[j]][[i]],na.rm =TRUE)
      desvpad <- sd(lista_das_listas_mpa_T[[j]][[i]],na.rm =TRUE)
      median <- median(lista_das_listas_mpa_T[[j]][[i]],na.rm =TRUE)
      n_origem <- n_orig[[j]]
      coeficiente <- alfa_h1[[i]]
      
      aux<-data.frame(media,
                      desvpad,
                      median,
                      n_origem,
                      coeficiente)
      df_resultados_alT<-rbind(df_resultados_alT,aux)
      
    }
   
  }
}


#  alt Q ------------------------------------------

for(j in 1:length(lista_das_listas_mpa_Q)){
  for( i in 1:length(lista_alfa_drift)){
 
    if(i==1 & j== 1){
      
      
      media <- mean(lista_das_listas_mpa_Q[[j]][[i]],na.rm =TRUE)
      desvpad <- sd(lista_das_listas_mpa_Q[[j]][[i]],na.rm =TRUE)
      median <- median(lista_das_listas_mpa_Q[[j]][[i]],na.rm =TRUE)
      n_origem <- n_orig[[j]]
      coeficiente <- alfa_h1[[i]]
      
      df_resultados_alQ<-data.frame(media,
                                    desvpad,
                                    median,
                                    n_origem,
                                    coeficiente)
      
    }else{
      
      aux<-NULL
      media <- mean(lista_das_listas_mpa_Q[[j]][[i]],na.rm =TRUE)
      desvpad <- sd(lista_das_listas_mpa_Q[[j]][[i]],na.rm =TRUE)
      median <- median(lista_das_listas_mpa_Q[[j]][[i]],na.rm =TRUE)
      n_origem <- n_orig[[j]]
      coeficiente <- alfa_h1[[i]]
      
      aux<-data.frame(media,
                      desvpad,
                      median,
                      n_origem,
                      coeficiente)
      df_resultados_alQ<-rbind(df_resultados_alQ,aux)
      
    }
  }
}

Vs    <-df_resultados_acf %>% dplyr::arrange(coeficiente, n_origem)
Vd    <-df_resultados_drift %>% dplyr::arrange(coeficiente, n_origem)
T2pca<-df_resultados_alT %>% dplyr::arrange(coeficiente, n_origem)
Q    <-df_resultados_alQ %>% dplyr::arrange(coeficiente, n_origem)
Vs 
Vd
T2pca
Q 
