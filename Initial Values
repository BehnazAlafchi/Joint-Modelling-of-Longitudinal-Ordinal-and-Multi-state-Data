

###################################################################################################################################################
################## Libraries ######################################################################################################################
###################################################################################################################################################
library("mixor")
library("mstate")
library("reshape2")
library("data.table")
###################################################################################################################################################

###################################################################################################################################################
################ Initial values ###################################################################################################################
###################################################################################################################################################
Initial<-function(data,associations,landa5points){
  init<-c()
  g<-c()
  beta<-c()
  sigma<-NULL
  alpha12<-c()
  alpha13<-c()
  alpha23<-c()
  
  datal<-melt(setDT(data), measure.vars=list(c(8,9,10),c(11,12,13)), 
              variable.name='var', value.name=c('t', 'y'))[,var:= paste0('measure',var)][order(id)]
  datal<-datal[,-2]
  datal<-cbind(datal, tx1=datal$t*datal$x1)
  
  l <- mixor(as.numeric(y) ~ t+x1+tx1, data = datal, which.random.slope = 1, id =id, link = "logit")
  lbeta<-c(l$Model[,1][c(3,2,4)])
  lsigma<-matrix(c(l$Model[,1][c(5,6,6,7)]),2,2)
  cm <- matrix(c(-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 1), ncol = 3)
  gama<-Contrasts(l, contrast.matrix = cm)
  lg<-gama[,1]
  
  tmat <- matrix(NA, 3, 3)
  tmat[1, 2:3] <- 1:2
  tmat[2, 3] <- 3
  dimnames(tmat) <- list(from = c("state0", "state1", "state2"), to = c("state0", "state1", "state2"))
  tmat
  covs <- c("x1")
  msbmt <- msprep(time = c(NA, "tir", "tid"), status = c(NA,"sir", "sid"), data = exam, trans = tmat,  keep = covs)
  msbmt <- expand.covs(msbmt, covs, append = TRUE, longnames = FALSE)
  head(msbmt)
  m <- coxph(Surv(Tstart, Tstop, status) ~ x1.1+x1.2+x1.3, data = msbmt, method = "breslow")
  mcoef<-m$coefficients
  malpha12=mcoef[1]
  malpha23=mcoef[3]
  malpha13=mcoef[2]
  
  g<-c(g,lg)
  beta<-c(beta, lbeta)
  sigma<-lsigma
  alpha12<-c(alpha12, malpha12)
  alpha13<-c(alpha13, malpha13)
  alpha23<-c(alpha23, malpha23)
  teta1=lsigma[1,1]
  teta2=lsigma[2,2]
  ro=((lsigma[1,2])/sqrt(lsigma[1,1]*lsigma[2,2]))
  param<-c(landa12=as.vector(landa5points),landa23=as.vector(landa5points),landa13=as.vector(landa5points),alpha12=alpha12
           ,alpha23=alpha23,alpha13=alpha13,eta12=associations[1]
           ,eta23=associations[2],eta13=associations[3]
           ,thre=g,beta=beta
           ,teta1=teta1,teta2=teta2
           ,ro=ro)
  return(param)
}

###################################################################################################################################################
