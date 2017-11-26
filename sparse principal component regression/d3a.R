####simulation####
library(mvtnorm)#use this package to generate normal dist
library(spcr)#use the packages to speed up this boring task
options(digits=3)#for simplicity
set.seed(666)
setwd("D:/study/DS2-1/stat760/final")
#case 3(a)
data3a = function(){
  MSE_1_50_0.1 = numeric(100)
  P_1_50_0.1 = numeric(100)
  N_1_50_0.1 = numeric(100)
  MSEA_1_50_0.1 = numeric(100)
  AP_1_50_0.1 = numeric(100)
  AN_1_50_0.1 = numeric(100)
  MSE_1_200_0.1 = numeric(100)
  P_1_200_0.1 = numeric(100)
  N_1_200_0.1 = numeric(100)
  MSEA_1_200_0.1 = numeric(100)
  AP_1_200_0.1 = numeric(100)
  AN_1_200_0.1 = numeric(100)
  MSE_10_50_0.1 = numeric(100)
  P_10_50_0.1 = numeric(100)
  N_10_50_0.1 = numeric(100)
  MSEA_10_50_0.1 = numeric(100)
  AP_10_50_0.1 = numeric(100)
  AN_10_50_0.1 = numeric(100)
  MSE_10_200_0.1 = numeric(100)
  P_10_200_0.1 = numeric(100)
  N_10_200_0.1 = numeric(100)
  MSEA_10_200_0.1 = numeric(100)
  AP_10_200_0.1 = numeric(100)
  AN_10_200_0.1 = numeric(100)
  MSE_1_50_1 = numeric(100)
  P_1_50_1 = numeric(100)
  N_1_50_1 = numeric(100)
  MSEA_1_50_1 = numeric(100)
  AP_1_50_1 = numeric(100)
  AN_1_50_1 = numeric(100)
  MSE_1_200_1 = numeric(100)
  P_1_200_1 = numeric(100)
  N_1_200_1 = numeric(100)
  MSEA_1_200_1 = numeric(100)
  AP_1_200_1 = numeric(100)
  AN_1_200_1 = numeric(100)
  MSE_10_50_1 = numeric(100)
  P_10_50_1 = numeric(100)
  N_10_50_1 = numeric(100)
  MSEA_10_50_1 = numeric(100)
  AP_10_50_1 = numeric(100)
  AN_10_50_1 = numeric(100)
  MSE_10_200_1 = numeric(100)
  P_10_200_1 = numeric(100)
  N_10_200_1 = numeric(100)
  MSEA_10_200_1 = numeric(100)
  AP_10_200_1 = numeric(100)
  AN_10_200_1 = numeric(100)
  for(t in 1:100){
    print(paste("this is",t,"simulation"))
    nonzero = c(T,F,T,T,F,T,T,F,T,rep(TRUE,6),rep(FALSE,15))
    error1 = rnorm(n=1250,mean = 0,sd=0.1)
    error2 = rnorm(n=1250,mean = 0,sd=1)
    sigma21=diag(rep(0,9))
    sigma31=diag(rep(0,6))
    for(i in 1:6)
    {
      for(j in 1:6)
      {
        sigma31[i,j]=0.9^abs(i-j)
      }
    }
    sigma3=cbind(rbind(sigma21,matrix(0,nrow=21,ncol=9)),rbind(matrix(0,nrow=9,ncol=6),sigma31,
                                                               matrix(0,nrow=15,ncol=6)),rbind(matrix(0,nrow=15,ncol=15),diag(15)))
    x3 = rmvnorm(n=1250,mean = rep(0,30),sigma = sigma3,method="chol")
    beta31=rep(1,6)
    beta3a=c(-1,0,1,1,0,-1,-1,0,1,beta31,rep(0,15))
    y3a1 = 4*x3%*%beta3a+error1
    y3a2 = 4*x3%*%beta3a+error2
    x = x3
    y1 = y3a1
    y2 = y3a2
    cv = cv.spcr(x=x[1:50,], y=y1[1:50], k=1)
    m = spcr(x=x[1:50,],y=y1[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    N_1_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[1:50,],y=y1[1:50],k=1,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y1[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    AN_1_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[51:250,], y=y1[51:250], k=1)
    m = spcr(x=x[51:250,],y=y1[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    N_1_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[51:250,],y=y1[51:250],k=1,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y1[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    AN_1_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[1:50,], y=y1[1:50], k=10)
    m = spcr(x=x[1:50,],y=y1[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    N_10_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[1:50,],y=y1[1:50],k=10,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y1[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    AN_10_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[51:250,], y=y1[51:250], k=10)
    m = spcr(x=x[51:250,],y=y1[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    N_10_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[51:250,],y=y1[51:250],k=10,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y1[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    AN_10_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[1:50,], y=y2[1:50], k=1)
    m = spcr(x=x[1:50,],y=y2[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    N_1_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[1:50,],y=y2[1:50],k=1,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y2[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    AN_1_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[51:250,], y=y2[51:250], k=1)
    m = spcr(x=x[51:250,],y=y2[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    N_1_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[51:250,],y=y2[51:250],k=1,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y2[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    AN_1_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[1:50,], y=y2[1:50], k=10)
    m = spcr(x=x[1:50,],y=y2[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    N_10_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[1:50,],y=y2[1:50],k=10,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y2[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    AN_10_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[51:250,], y=y2[51:250], k=10)
    m = spcr(x=x[51:250,],y=y2[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    N_10_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
    cv = cv.spcr(x=x[51:250,],y=y2[51:250],k=10,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y2[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/12
    AN_10_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/18
  }
  result = rbind(MSE_1_50_0.1,P_1_50_0.1,N_1_50_0.1,AP_1_50_0.1,AN_1_50_0.1,MSEA_1_50_0.1,P_1_50_0.1,N_1_50_0.1,
                 AP_1_50_0.1,AN_1_50_0.1,MSE_1_200_0.1,P_1_200_0.1,N_1_200_0.1,AP_1_200_0.1,AN_1_200_0.1,
                 MSEA_1_200_0.1,P_1_200_0.1,N_1_200_0.1,AP_1_200_0.1,AN_1_200_0.1,MSE_10_50_0.1,P_10_50_0.1,N_10_50_0.1,
                 AP_10_50_0.1,AN_10_50_0.1,MSEA_10_50_0.1,P_10_50_0.1,N_10_50_0.1,AP_10_50_0.1,AN_10_50_0.1,
                 MSE_10_200_0.1,P_10_200_0.1,N_10_200_0.1,AP_10_200_0.1,AN_10_200_0.1,MSEA_10_200_0.1,P_10_200_0.1,
                 N_10_200_0.1,AP_10_200_0.1,AN_10_200_0.1,MSE_1_50_1,P_1_50_1,N_1_50_1,AP_1_50_1,AN_1_50_1,MSEA_1_50_1,
                 P_1_50_1,N_1_50_1,AP_1_50_1,AN_1_50_1,MSE_1_200_1,P_1_200_1,N_1_200_1,AP_1_200_1,AN_1_200_1,MSEA_1_200_1,
                 P_1_200_1,N_1_200_1,AP_1_200_1,AN_1_200_1,MSE_10_50_1,P_10_50_1,N_10_50_1,AP_10_50_1,AN_10_50_1,MSEA_10_50_1,
                 P_10_50_1,N_10_50_1,AP_10_50_1,AN_10_50_1,MSE_10_200_1,P_10_200_1,N_10_200_1,AP_10_200_1,AN_10_200_1,MSEA_10_200_1,
                 P_10_200_1,N_10_200_1,AP_10_200_1,AN_10_200_1)
  
  MEAN = apply(result,1,mean)
  SD = apply(result,1,sd)
  ans = list(mean = MEAN,sd = SD)
  return(ans)
}
ans = data3a()
haha = cbind(mean = ans$mean,sd = ans$sd)
write.table(haha,"data3a.csv")