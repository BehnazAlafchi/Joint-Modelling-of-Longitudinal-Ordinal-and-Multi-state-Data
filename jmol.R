jmol<-function(data,associations,landa5points){
  Result<-list()
  like<-function(par,data){
    #death
    landa31 <- exp(par[1])
    landa32 <- exp(par[2])
    landa33 <- exp(par[3])
    landa34 <- exp(par[4])
    landa35 <- exp(par[5])
    
    #recurrent
    landa11 <- exp(par[6])
    landa12 <- exp(par[7])
    landa13 <- exp(par[8])
    landa14 <- exp(par[9])
    landa15 <- exp(par[10])
    
    #gap
    landa21 <- exp(par[11])
    landa22 <- exp(par[12])
    landa23 <- exp(par[13])
    landa24 <- exp(par[14])
    landa25 <- exp(par[15])
    
    
    alpha12=par[16]
    alpha23=par[17]
    alpha13=par[18]
    eta12=par[19]
    eta23=par[20]
    eta13=par[21]
    g=par[22:24]
    beta=par[25:27]
    teta1=exp(par[28])
    teta2=exp(par[29])
    r=par[30]
    
    tid=data$tid
    tir=data$tir
    sid=data$sid
    sir=data$sir
    sig-data$sig
    t_1=data$t_1
    t_2=data$t_2
    t_3=data$t_3
    y_1=data$y_1
    y_2=data$y_2
    y_3=data$y_3
    x1=data$x1
    
    tid=tid+0.01
    tig=tid-tir
    
    
    y_1=as.numeric(y_1)
    y_2=as.numeric(y_2)
    y_3=as.numeric(y_3)
    
    t=matrix(c(t_1,t_2,t_3),nrow(data),3)
    y=matrix(c(y_1,y_2,y_3),nrow(data),3)
    x=matrix(c(x1),nrow(data),2)
    
    
    myGrid1 <- createNIGrid(dim=2, type="GHe", level=5)
    
    nn=3
    lqq<-c()
    lvv<-c()
    lv<-c()
    deter=qq=v=vv=0
    
    baseHM <- baseHR<- baseHG <- cumbaseM <- cumbaseR<- cumbaseG <- rep(1, nrow(data))
    tQD <- as.vector(quantile(tid[sid == 1], probs = c(0.2, 0.4, 0.6, 0.8, 1)))
    tQR <- as.vector(quantile(tir[sir == 1], probs = c(0.2, 0.4, 0.6, 0.8, 1)))
    tQG <- as.vector(quantile(tig[sig == 1], probs = c(0.2, 0.4, 0.6, 0.8, 1)))
    
    
    for(i in 1:nrow(data)){
      
      
      ####################  #######################  
      
      
      #longitudinal
      #m=the number of repeated measurements
      m=3
      #c=the number of level
      c=4
      long<-function(g,beta,b){
        a<-NULL
        a0<-1
        for(j in 1:m){
          if(y[i,j]==1) a<-(1/(1+exp(-(g[y[i,j]]-(sum(beta*c(x[i,],t[i,j],(x[i,]*t[i,j])))+b[,1]+(b[,2]*t[i,j]))))))-(0)
          if(y[i,j]>1 & y[i,j]<c) a<-(1/(1+exp(-(g[y[i,j]]-(sum(beta*c(x[i,],t[i,j],(x[i,]*t[i,j])))+b[,1]+(b[,2]*t[i,j]))))))-(1/(1+exp(-(g[y[i,j]-1]-(sum(beta*c(x[i,],t[i,j],(x[i,]*t[i,j])))+b[,1]+(b[,2]*t[i,j]))))))
          if(y[i,j]==c) a<-(1)-(1/(1+exp(-(g[y[i,j]-1]-(sum(beta*c(x[i,],t[i,j],(x[i,]*t[i,j])))+b[,1]+(b[,2]*t[i,j]))))))
          a0<-a0*a
        }
        return(a0)
      }
      
      #Random effects
      rand <- function(b,teta1,teta2,r) {
        exp(-((1/(2*(1-r)))*(((b[,1]^2)/teta1)-((2*r*b[,1]*b[,2])/(sqrt(teta1*teta2)))+((b[,2]^2)/teta2))))
      }
      
      #multi-state (including quadrature)
      #1
      if (sir[i] == 1 & sid[i] == 1){
        
        deter=sum(alpha12*x[i,])+sum(alpha23*x[i,])-0.5*log(2*pi)-0.5*log(teta1*teta2*(1-r^2))
        v=deter
        multi<-function(eta12,eta23,b,landa11,landa12,landa13,landa14,landa15,landa21, landa22,landa23,landa24,landa25,alpha12,alpha23){
          
          if ( tir[i] > 0 & tir[i] <= tQR[1]) {
            baseHR[i] <- landa11*exp(eta12*b[2]*tir[i])
            cumbaseR[i] <- ((exp(eta12*b[1]))/(eta12*b[2]))*(landa11*(exp(eta12*b[2]*tir[i])-1))
          }
          if (tir[i] > tQR[1] & tir[i] <= tQR[2]) {
            baseHR[i] <- landa12*exp(eta12*b[2]*tir[i])
            cumbaseR[i] <- ((exp(eta12*b[1]))/(eta12*b[2]))* ((landa11*((exp(eta12*b[2]*tQR[1]))-1))+ (landa12 * (exp(eta12*b[2]*tir[i])-exp(eta12*b[2]*tQR[1]))))
          }
          if (tir[i] > tQR[2] & tir[i] <= tQR[3]) {
            baseHR[i] <- landa13*exp(eta12*b[2]*tir[i])
            cumbaseR[i] <- ((exp(eta12*b[1]))/(eta12*b[2]))*((landa11*((exp(eta12*b[2]*tQR[1]))-1))+ (landa12 * (exp(eta12*b[2]*tQR[2])-exp(eta12*b[2]*tQR[1])))
                                                             + (landa13 * (exp(eta12*b[2]*tir[i])-exp(eta12*b[2]*tQR[2]))) )
          }
          if (tir[i] > tQR[3] & tir[i] <=  tQR[4]) {
            baseHR[i] <- landa14*exp(eta12*b[2]*tir[i])
            cumbaseR[i] <- ((exp(eta12*b[1]))/(eta12*b[2]))*((landa11*((exp(eta12*b[2]*tQR[1]))-1))+ (landa12 * (exp(eta12*b[2]*tQR[2])-exp(eta12*b[2]*tQR[1])))
                                                             + (landa13 * (exp(eta12*b[2]*tQR[3])-exp(eta12*b[2]*tQR[2])))  + (landa14 * (exp(eta12*b[2]*tir[i])-exp(eta12*b[2]*tQR[3]))) )
          }
          if (tir[i] > tQR[4]) {
            baseHR[i] <- landa15*exp(eta12*b[2]*tir[i])
            cumbaseR[i] <- ((exp(eta12*b[1]))/(eta12*b[2]))*((landa11*((exp(eta12*b[2]*tQR[1]))-1))+ (landa12 * (exp(eta12*b[2]*tQR[2])-exp(eta12*b[2]*tQR[1])))
                                                             + (landa13 * (exp(eta12*b[2]*tQR[3])-exp(eta12*b[2]*tQR[2])))  + (landa14 * (exp(eta12*b[2]*tQR[4])-exp(eta12*b[2]*tQR[3]))) 
                                                             + (landa15 * (exp(eta12*b[2]*tir[i])-exp(eta12*b[2]*tQR[4]))) )
          }
          if ( tig[i] > 0 & tig[i] <= tQG[1]) {
            baseHG[i] <- landa21*exp(eta23*b[2]*tig[i])
            cumbaseG[i] <- ((exp(eta23*b[1]))/(eta23*b[2]))*(landa21*(exp(eta23*b[2]*tig[i])-1))
          }
          if (tig[i] > tQG[1] & tig[i] <= tQG[2]) {
            baseHG[i] <- landa22*exp(eta23*b[2]*tig[i])
            cumbaseG[i] <- ((exp(eta23*b[1]))/(eta23*b[2]))* ((landa21*((exp(eta23*b[2]*tQG[1]))-1))+ (landa22 * (exp(eta23*b[2]*tig[i])-exp(eta23*b[2]*tQG[1]))))
          }
          if (tig[i] > tQG[2] & tig[i] <= tQG[3]) {
            baseHG[i] <- landa23*exp(eta23*b[2]*tig[i])
            cumbaseG[i] <- ((exp(eta23*b[1]))/(eta23*b[2]))*((landa21*((exp(eta23*b[2]*tQG[1]))-1))+ (landa22 * (exp(eta23*b[2]*tQG[2])-exp(eta23*b[2]*tQG[1])))
                                                             + (landa23 * (exp(eta23*b[2]*tig[i])-exp(eta23*b[2]*tQG[2]))) )
          }
          if (tig[i] > tQG[3] & tig[i] <=  tQG[4]) {
            baseHG[i] <- landa24*exp(eta23*b[2]*tig[i])
            cumbaseG[i] <- ((exp(eta23*b[1]))/(eta23*b[2]))*((landa21*((exp(eta23*b[2]*tQG[1]))-1))+ (landa22 * (exp(eta23*b[2]*tQG[2])-exp(eta23*b[2]*tQG[1])))
                                                             + (landa23 * (exp(eta23*b[2]*tQG[3])-exp(eta23*b[2]*tQG[2])))  + (landa24 * (exp(eta23*b[2]*tig[i])-exp(eta23*b[2]*tQG[3]))) )
          }
          if (tig[i] > tQG[4]) {
            baseHG[i] <- landa25*exp(eta23*b[2]*tig[i])
            cumbaseG[i] <- ((exp(eta23*b[1]))/(eta23*b[2]))*((landa21*((exp(eta23*b[2]*tQG[1]))-1))+ (landa22 * (exp(eta23*b[2]*tQG[2])-exp(eta23*b[2]*tQG[1])))
                                                             + (landa23 * (exp(eta23*b[2]*tQG[3])-exp(eta23*b[2]*tQG[2])))  + (landa24 * (exp(eta23*b[2]*tQG[4])-exp(eta23*b[2]*tQG[3]))) 
                                                             + (landa25 * (exp(eta23*b[2]*tig[i])-exp(eta23*b[2]*tQG[4]))) )
          }
          
          ((exp(-(exp(sum(alpha12*x[i,])))* cumbaseR[i]))*(exp(eta12*b[1]))*baseHR[i])*
            (exp(-(exp(sum(alpha23*x[i,])))* cumbaseG[i]))*(exp(eta23*b[1]))*baseHG[i]
          
        }
        myfun<-function(b,g,beta,eta12,eta23,landa11,landa12,landa13,landa14,landa15,landa21, landa22,landa23,landa24,landa25
                        ,alpha12,alpha23,teta1,teta2,r){
          long(g,beta,b)*multi(eta12,eta23,b,landa11,landa12,landa13,landa14,landa15,landa21, landa22,landa23,landa24,landa25,
                               alpha12,alpha23)*rand (b,teta1,teta2,r)
        }
        qq<- quadrature(f = myfun,g=g,beta=beta,eta12=eta12,eta23=eta23,landa11=landa11,landa12=landa12,landa13=landa13, landa14=landa14, 
                        landa15=landa15,landa21=landa21,landa22=landa22,landa23=landa23,landa24=landa24,landa25=landa25,alpha12=alpha12,
                        alpha23=alpha23,teta1=teta1,teta2=teta2,r=r,grid = myGrid1)
        v<-v+log(qq)
      }
      #2
      if (sir[i] == 1 & sid[i] == 0){
        
        deter=sum(alpha12*x[i,])-0.5*log(2*pi)-0.5*log(teta1*teta2*(1-r^2))
        v=deter
        multi<-function(eta12,eta23,b,landa11,landa12,landa13,landa14,landa15,landa21, landa22,landa23,landa24,landa25 ,alpha12,alpha23){
          
          
          if ( tir[i] > 0 & tir[i] <= tQR[1]) {
            baseHR[i] <- landa11*exp(eta12*b[2]*tir[i])
            cumbaseR[i] <- ((exp(eta12*b[1]))/(eta12*b[2]))*(landa11*(exp(eta12*b[2]*tir[i])-1))
          }
          if (tir[i] > tQR[1] & tir[i] <= tQR[2]) {
            baseHR[i] <- landa12*exp(eta12*b[2]*tir[i])
            cumbaseR[i] <- ((exp(eta12*b[1]))/(eta12*b[2]))* ((landa11*((exp(eta12*b[2]*tQR[1]))-1))+ (landa12 * (exp(eta12*b[2]*tir[i])-exp(eta12*b[2]*tQR[1]))))
          }
          if (tir[i] > tQR[2] & tir[i] <= tQR[3]) {
            baseHR[i] <- landa13*exp(eta12*b[2]*tir[i])
            cumbaseR[i] <- ((exp(eta12*b[1]))/(eta12*b[2]))*((landa11*((exp(eta12*b[2]*tQR[1]))-1))+ (landa12 * (exp(eta12*b[2]*tQR[2])-exp(eta12*b[2]*tQR[1])))
                                                             + (landa13 * (exp(eta12*b[2]*tir[i])-exp(eta12*b[2]*tQR[2]))) )
          }
          if (tir[i] > tQR[3] & tir[i] <=  tQR[4]) {
            baseHR[i] <- landa14*exp(eta12*b[2]*tir[i])
            cumbaseR[i] <- ((exp(eta12*b[1]))/(eta12*b[2]))*((landa11*((exp(eta12*b[2]*tQR[1]))-1))+ (landa12 * (exp(eta12*b[2]*tQR[2])-exp(eta12*b[2]*tQR[1])))
                                                             + (landa13 * (exp(eta12*b[2]*tQR[3])-exp(eta12*b[2]*tQR[2])))  + (landa14 * (exp(eta12*b[2]*tir[i])-exp(eta12*b[2]*tQR[3]))) )
          }
          if (tir[i] > tQR[4]) {
            baseHR[i] <- landa15*exp(eta12*b[2]*tir[i])
            cumbaseR[i] <- ((exp(eta12*b[1]))/(eta12*b[2]))*((landa11*((exp(eta12*b[2]*tQR[1]))-1))+ (landa12 * (exp(eta12*b[2]*tQR[2])-exp(eta12*b[2]*tQR[1])))
                                                             + (landa13 * (exp(eta12*b[2]*tQR[3])-exp(eta12*b[2]*tQR[2])))  + (landa14 * (exp(eta12*b[2]*tQR[4])-exp(eta12*b[2]*tQR[3]))) 
                                                             + (landa15 * (exp(eta12*b[2]*tir[i])-exp(eta12*b[2]*tQR[4]))) )
          }
          if ( tig[i] > 0 & tig[i] <= tQG[1]) {
            baseHG[i] <- landa21*exp(eta23*b[2]*tig[i])
            cumbaseG[i] <- ((exp(eta23*b[1]))/(eta23*b[2]))*(landa21*(exp(eta23*b[2]*tig[i])-1))
          }
          if (tig[i] > tQG[1] & tig[i] <= tQG[2]) {
            baseHG[i] <- landa22*exp(eta23*b[2]*tig[i])
            cumbaseG[i] <- ((exp(eta23*b[1]))/(eta23*b[2]))* ((landa21*((exp(eta23*b[2]*tQG[1]))-1))+ (landa22 * (exp(eta23*b[2]*tig[i])-exp(eta23*b[2]*tQG[1]))))
          }
          if (tig[i] > tQG[2] & tig[i] <= tQG[3]) {
            baseHG[i] <- landa23*exp(eta23*b[2]*tig[i])
            cumbaseG[i] <- ((exp(eta23*b[1]))/(eta23*b[2]))*((landa21*((exp(eta23*b[2]*tQG[1]))-1))+ (landa22 * (exp(eta23*b[2]*tQG[2])-exp(eta23*b[2]*tQG[1])))
                                                             + (landa23 * (exp(eta23*b[2]*tig[i])-exp(eta23*b[2]*tQG[2]))) )
          }
          if (tig[i] > tQG[3] & tig[i] <=  tQG[4]) {
            baseHG[i] <- landa24*exp(eta23*b[2]*tig[i])
            cumbaseG[i] <- ((exp(eta23*b[1]))/(eta23*b[2]))*((landa21*((exp(eta23*b[2]*tQG[1]))-1))+ (landa22 * (exp(eta23*b[2]*tQG[2])-exp(eta23*b[2]*tQG[1])))
                                                             + (landa23 * (exp(eta23*b[2]*tQG[3])-exp(eta23*b[2]*tQG[2])))  + (landa24 * (exp(eta23*b[2]*tig[i])-exp(eta23*b[2]*tQG[3]))) )
          }
          if (tig[i] > tQG[4]) {
            baseHG[i] <- landa25*exp(eta23*b[2]*tig[i])
            cumbaseG[i] <- ((exp(eta23*b[1]))/(eta23*b[2]))*((landa21*((exp(eta23*b[2]*tQG[1]))-1))+ (landa22 * (exp(eta23*b[2]*tQG[2])-exp(eta23*b[2]*tQG[1])))
                                                             + (landa23 * (exp(eta23*b[2]*tQG[3])-exp(eta23*b[2]*tQG[2])))  + (landa24 * (exp(eta23*b[2]*tQG[4])-exp(eta23*b[2]*tQG[3]))) 
                                                             + (landa25 * (exp(eta23*b[2]*tig[i])-exp(eta23*b[2]*tQG[4]))) )
          }
          
          
          ((exp(-(exp(sum(alpha12*x[i,])))* cumbaseR[i]))*(exp(eta12*b[1]))*baseHR[i])*
            (exp(-(exp(sum(alpha23*x[i,])))* cumbaseG[i]))
          
        }
        myfun<-function(b,g,beta,eta12,eta23,landa11,landa12,landa13,landa14,landa15,landa21, landa22,landa23,landa24,landa25,
                        alpha12,alpha23,teta1,teta2,r){
          long(g,beta,b)*multi(eta12,eta23,b,landa11,landa12,landa13,landa14,landa15,landa21, landa22,landa23,landa24,landa25,
                               alpha12,alpha23)*rand (b,teta1,teta2,r)
        }
        qq<- quadrature(f = myfun,g=g,beta=beta,eta12=eta12,eta23=eta23,landa11=landa11,landa12=landa12,landa13=landa13, landa14=landa14, 
                        landa15=landa15, landa21=landa21,landa22=landa22,landa23=landa23,landa24=landa24,landa25=landa25,alpha12=alpha12,alpha23=alpha23,teta1=teta1,teta2=teta2,r=r,grid = myGrid1)
        v<-v+log(qq)
      }
      #3
      if (sir[i] == 0 & sid[i] == 1){
        
        deter=sum(alpha13*x[i,])-0.5*log(2*pi)-0.5*log(teta1*teta2*(1-r^2))
        v=deter
        multi<-function(eta13,b,landa31,landa32, landa33, landa34, landa35,alpha13){
          
          if ( tid[i] > 0 & tid[i] <= tQD[1]) {
            baseHM[i] <- landa31*exp(eta13*b[2]*tid[i])
            cumbaseM[i] <- ((exp(eta13*b[1]))/(eta13*b[2]))*(landa31*(exp(eta13*b[2]*tid[i])-1))
          }
          if (tid[i] > tQD[1] & tid[i] <= tQD[2]) {
            baseHM[i] <- landa32*exp(eta13*b[2]*tid[i])
            cumbaseM[i] <- ((exp(eta13*b[1]))/(eta13*b[2]))* ((landa31*((exp(eta13*b[2]*tQD[1]))-1))+ (landa32 * (exp(eta13*b[2]*tid[i])-exp(eta13*b[2]*tQD[1]))))
          }
          if (tid[i] > tQD[2] & tid[i] <= tQD[3]) {
            baseHM[i] <- landa33*exp(eta13*b[2]*tid[i])
            cumbaseM[i] <- ((exp(eta13*b[1]))/(eta13*b[2]))*((landa31*((exp(eta13*b[2]*tQD[1]))-1))+ (landa32 * (exp(eta13*b[2]*tQD[2])-exp(eta13*b[2]*tQD[1])))
                                                             + (landa33 * (exp(eta13*b[2]*tid[i])-exp(eta13*b[2]*tQD[2]))) )
          }
          if (tid[i] > tQD[3] & tid[i] <=  tQD[4]) {
            baseHM[i] <- landa34*exp(eta13*b[2]*tid[i])
            cumbaseM[i] <- ((exp(eta13*b[1]))/(eta13*b[2]))*((landa31*((exp(eta13*b[2]*tQD[1]))-1))+ (landa32 * (exp(eta13*b[2]*tQD[2])-exp(eta13*b[2]*tQD[1])))
                                                             + (landa33 * (exp(eta13*b[2]*tQD[3])-exp(eta13*b[2]*tQD[2])))  + (landa34 * (exp(eta13*b[2]*tid[i])-exp(eta13*b[2]*tQD[3]))) )
          }
          if (tid[i] > tQD[4]) {
            baseHM[i] <- landa35*exp(eta13*b[2]*tid[i])
            cumbaseM[i] <- ((exp(eta13*b[1])) /(eta13*b[2]))*((landa31*((exp(eta13*b[2]*tQD[1]))-1))+ (landa32 * (exp(eta13*b[2]*tQD[2])-exp(eta13*b[2]*tQD[1])))
                                                              + (landa33 * (exp(eta13*b[2]*tQD[3])-exp(eta13*b[2]*tQD[2])))  + (landa34 * (exp(eta13*b[2]*tQD[4])-exp(eta13*b[2]*tQD[3]))) 
                                                              + (landa35 * (exp(eta13*b[2]*tid[i])-exp(eta13*b[2]*tQD[4]))) )
          }
          (exp(-(exp(sum(alpha13*x[i,])))* cumbaseM[i]))*(exp(eta13*b[1]))*baseHM[i]
        }
        myfun<-function(b,g,beta,eta13,landa31,landa32, landa33, landa34, landa35,alpha13,teta1,teta2,r){
          long(g,beta,b)*multi(eta13,b,landa31,landa32, landa33, landa34, landa35,alpha13)*rand (b,teta1,teta2,r)
        }
        qq<- quadrature(f = myfun,g=g,beta=beta,eta13=eta13,landa31=landa31,landa32=landa32, landa33=landa33, 
                        landa34=landa34, landa35=landa35,alpha13=alpha13,teta1=teta1,teta2=teta2,r=r,grid = myGrid1)
        v<-v+log(qq)
      }
      #4
      if (sir[i] == 0 & sid[i] == 0){
        
        deter=-0.5*log(2*pi)-0.5*log(teta1*teta2*(1-r^2))
        v=deter
        multi<-function(eta13,b,landa31,landa32, landa33, landa34, landa35,alpha13){
          
          if ( tid[i] > 0 & tid[i] <= tQD[1]) {
            baseHM[i] <- landa31*exp(eta13*b[2]*tid[i])
            cumbaseM[i] <- ((exp(eta13*b[1]))/(eta13*b[2]))*(landa31*(exp(eta13*b[2]*tid[i])-1))
          }
          if (tid[i] > tQD[1] & tid[i] <= tQD[2]) {
            baseHM[i] <- landa32*exp(eta13*b[2]*tid[i])
            cumbaseM[i] <- ((exp(eta13*b[1]))/(eta13*b[2]))* ((landa31*((exp(eta13*b[2]*tQD[1]))-1))+ (landa32 * (exp(eta13*b[2]*tid[i])-exp(eta13*b[2]*tQD[1]))))
          }
          if (tid[i] > tQD[2] & tid[i] <= tQD[3]) {
            baseHM[i] <- landa33*exp(eta13*b[2]*tid[i])
            cumbaseM[i] <- ((exp(eta13*b[1]))/(eta13*b[2]))*((landa31*((exp(eta13*b[2]*tQD[1]))-1))+ (landa32 * (exp(eta13*b[2]*tQD[2])-exp(eta13*b[2]*tQD[1])))
                                                             + (landa33 * (exp(eta13*b[2]*tid[i])-exp(eta13*b[2]*tQD[2]))) )
          }
          if (tid[i] > tQD[3] & tid[i] <=  tQD[4]) {
            baseHM[i] <- landa34*exp(eta13*b[2]*tid[i])
            cumbaseM[i] <- ((exp(eta13*b[1]))/(eta13*b[2]))*((landa31*((exp(eta13*b[2]*tQD[1]))-1))+ (landa32 * (exp(eta13*b[2]*tQD[2])-exp(eta13*b[2]*tQD[1])))
                                                             + (landa33 * (exp(eta13*b[2]*tQD[3])-exp(eta13*b[2]*tQD[2])))  + (landa34 * (exp(eta13*b[2]*tid[i])-exp(eta13*b[2]*tQD[3]))) )
          }
          if (tid[i] > tQD[4]) {
            baseHM[i] <- landa35*exp(eta13*b[2]*tid[i])
            cumbaseM[i] <- ((exp(eta13*b[1])) /(eta13*b[2]))*((landa31*((exp(eta13*b[2]*tQD[1]))-1))+ (landa32 * (exp(eta13*b[2]*tQD[2])-exp(eta13*b[2]*tQD[1])))
                                                              + (landa33 * (exp(eta13*b[2]*tQD[3])-exp(eta13*b[2]*tQD[2])))  + (landa34 * (exp(eta13*b[2]*tQD[4])-exp(eta13*b[2]*tQD[3]))) 
                                                              + (landa35 * (exp(eta13*b[2]*tid[i])-exp(eta13*b[2]*tQD[4]))) )
          }
          exp(-(exp(sum(alpha13*x[i,])))* cumbaseM[i])
          
        }
        
        myfun<-function(b,g,beta,eta13,landa31,landa32, landa33, landa34, landa35,alpha13,teta1,teta2,r){
          long(g,beta,b)*multi(eta13,b,landa31,landa32, landa33, landa34, landa35,alpha13)*rand (b,teta1,teta2,r)
        }
        qq<- quadrature(f = myfun,g=g,beta=beta,eta13=eta13,landa31=landa31,landa32=landa32, landa33=landa33, 
                        landa34=landa34, landa35=landa35,alpha13=alpha13,teta1=teta1,teta2=teta2,r=r,grid = myGrid1)
        v<-v+log(qq)
      }
      vv<-log(qq)
      lqq<-c(lqq,qq)
      lvv<-c(lvv,vv)
      lv<-c(lv,v)
    }
    return(total.likelihood=sum(lv))
  } 
  Inparam=Initial(data,associations,landa5points)

  likelihood<-like(Inparam,data=data)
  ml<-optim(Inparam, fn=like,data=data)
  ML<-ml$par[-c(1:15)]
  
  Result<-list(Coefficient=ML,likelihood_value=likelihood)
  return(Result)
}
