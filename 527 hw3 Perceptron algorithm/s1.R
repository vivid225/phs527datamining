s1 <- function(id){
  library(gridExtra)
  library(rjags)
  library(msm)
  library(VGAM) 
  library(robustbase)
  library(stats)
  library(ggplot2)
  library(betareg)
  library(boot)
  library(reshape2)
  library(xtable)
  library(gamlss)
  library(mgcv)
  library(drc)
  library(betaboost)
  
  set.seed(12)
  source("Betareg_functions_notruncate.R")
  source("MDPDE.r")
  source("SMLE.r")
  
  ## True dose-response curve -----
  
  d1=0.1 ;  dose.k=7 ;  multiplier=2 ;
  #generate 7 doses by d_k=d_1*(multiplier)^(k-1), 
  #where d_1 is arbitrary and multiplier chosen as 2 here
  
  dose1=d1
  dose.set<-NULL
  for (i in 1:(dose.k-1)) {
    dose.set[i]=dose1*multiplier^i
  }
  dose.set<-c(dose1, dose.set)
  
  
  ###given true dose-response curve  f(x) = inv.logit(b_0+b_1*log(x)) 
  ##let's have k doses, want f(d1)=.01, f(dk)=.99
  y1<-logit(0.01) #b_0+b_1*log(dose1)
  y2<-logit(0.99) #b_0+b_1*log(dose7)
  b_1=(y1-y2)/(log(dose.set[1])-log(dose.set[dose.k]))
  b_0=y1-b_1*log(dose.set[1]) 
  truef<-function(x){ 
    prop = inv.logit(b_0+b_1*log(x))
    return(prop)
  }
  
  ## simulate data Y -----
  n.sim =10
  phi = 10
  # id=1
  dose.set <- c(0.04, 0.08, 0.14, 0.27, 0.52, 0.99, 1.88, 3.58, 6.79, 12.91)
  # prop <- truef(dose.set)
  # propnew <- c(prop[1:6], sort(1-prop[1:4]))
  # exp((logit(propnew)-b_0)/b_1)
  
  ### 10 doses     
  Yy<-matrix(rep(0, length(dose.set)*n.sim*10), nrow=length(dose.set), ncol=n.sim*10)
  
  for (s in 1:length(dose.set)) {
    mu=truef(dose.set[s])
    Yy[s, ]=rbeta(n.sim, mu*phi, (1-mu)*phi)
  }
  
  Yy[Yy <= 0] = 1e-9
  Yy[Yy >= 1] = 1-1e-9
  
  # plot(Yy[,1],ylim=c(0,1));Yy[,1]
  
  logTTDm<-(logit(0.5)-b_0)/b_1
  TTDm<-exp(logTTDm)
  TTDm
  
  logTTD5<-(logit(0.05)-b_0)/b_1
  TTD5<-exp(logTTD5)
  TTD5
  
  logTTD95<-(logit(0.95)-b_0)/b_1
  TTD95<-exp(logTTD95)
  TTD95
  
  
  
  ## Running for results -----
  res_brm <- matrix(NA, n.sim, 15)
  res_brmp_lasso <- matrix(NA, n.sim, 5)
  res_brmp_ridge <- matrix(NA, n.sim, 5)
  res_brmq <- matrix(NA, n.sim, 15)
  res_mdpde <- matrix(NA, n.sim, 15)
  res_smle <- matrix(NA, n.sim, 15)
  res_t <- matrix(NA, n.sim, 15)
  res_drc <- matrix(NA, n.sim, 15)
  res_nls_gn <- matrix(NA, n.sim, 15)
  res_nls_gb <- matrix(NA, n.sim, 15)
  res_nls_port <- matrix(NA, n.sim, 15)
  res_betaboost <- matrix(NA, n.sim, 5)
  
  j = 1
  for(i in (1000*(id-1)+1):(1000*id)) {
    fit <- simuf_trun(Yy=Yy,i=i,philink="identity", vary_precision=FALSE,
                      BRM=TRUE,BRMP_LASSO=FALSE,BRMP_Ridge=TRUE,BRMQ=TRUE,
                      MDPDE=TRUE,SMLE=TRUE,BBRM_T=TRUE,
                      BBRM_G=FALSE,BBRM_C=FALSE, drc=TRUE,
                      nls_gn=TRUE,nls_gb=TRUE,nls_port=TRUE,
                      betaboost=TRUE)
    res_brm[j,] = fit$brm
    # res_brmp_lasso[i,] = fit$brmp_lasso
    res_brmp_ridge[j,] = fit$brmp_ridge
    res_brmq[j,] = fit$brmq
    res_mdpde[j,] = fit$mdpde
    res_smle[j,] = fit$smle
    res_t[j,] = fit$tprior
    res_drc[j,] = fit$drc
    res_nls_gn[j,] = fit$nls_gn
    res_nls_gb[j,] = fit$nls_gb
    res_nls_port[j,] = fit$nls_port
    res_betaboost[j,] = fit$betaboost
    j = j+1
  }
  
  restotal_identity <- list(brm = res_brm,
                            # brmp_lasso = res_brmp_lasso,
                            brmp_ridge = res_brmp_ridge,
                            brmq = res_brmq,
                            mdpde = res_mdpde,
                            smle = res_smle,
                            t = res_t,
                            drc = res_drc,
                            nls_gn = res_nls_gn,
                            nls_gb = res_nls_gb,
                            nls_port = res_nls_port,
                            betaboost = res_betaboost)
  
  rds <- paste("scenario1_",id,".rds",sep = "")
  saveRDS(restotal_identity, rds)
  
  # res_brm <- matrix(NA, n.sim, 15)
  # res_brmp_lasso <- matrix(NA, n.sim, 5)
  # res_brmp_ridge <- matrix(NA, n.sim, 5)
  # res_brmq <- matrix(NA, n.sim, 9)
  # res_mdpde <- matrix(NA, n.sim, 15)
  # res_smle <- matrix(NA, n.sim, 15)
  # res_t <- matrix(NA, n.sim, 15)
  # 
  # for(i in 1:n.sim) { 
  #   fit <- simuf_trun(i=i,philink="log", vary_precision=FALSE,
  #                     BRM=TRUE,BRMP=TRUE,BRMQ=FALSE,
  #                     MDPDE=TRUE,SMLE=TRUE,BBRM_T=TRUE,
  #                     BBRM_G=FALSE,BBRM_C=FALSE)
  #   res_brm[i,] = fit$brm
  #   res_brmp[i,] = fit$brmp
  #   res_mdpde[i,] = fit$mdpde
  #   res_smle[i,] = fit$smle
  #   res_t[i,] = fit$tprior
  # }   
  # 
  # restotal_log <- list(brm = res_brm,
  #                      brmp_lasso = res_brmp_lasso,
  #                      brmp_ridge = res_brmp_ridge,
  #                      brmq = res_brmq,
  #                      mdpde = res_mdpde,
  #                      smle = res_smle,
  #                      t = res_t)
  
  
  
  # saveRDS(restotal_log, "scenario11_11.rds")
  
  
  
}