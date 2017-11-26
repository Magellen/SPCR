####simulation####
library(mvtnorm)#use this package to generate normal dist
library(spcr)#use the packages to speed up this boring task
options(digits=3)#for simplicity
set.seed(666)
setwd("D:/study/DS2-1/stat760/final")
#case 1(b)
data1b = function(){
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
    nonzero = c(TRUE,TRUE,rep(FALSE,8))
    error1 = rnorm(n=1250,mean = 0,sd=0.1)
    error2 = rnorm(n=1250,mean = 0,sd=1)
    sigma1b=diag(rep(1,10))
    sigma1b[2,2]=9
    x1b = rmvnorm(n=1250,mean = rep(0,10),sigma = sigma1b,method="chol")
    beta1b = c(8,1,rep(0,8))
    y1b1 = x1b%*%beta1b+error1
    y1b2 = x1b%*%beta1b+error2
    x = x1b
    y1 = y1b1
    y2 = y1b2
    cv = cv.spcr(x=x[1:50,], y=y1[1:50], k=1)
    m = spcr(x=x[1:50,],y=y1[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    N_1_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[1:50,],y=y1[1:50],k=1,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y1[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    AN_1_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[51:250,], y=y1[51:250], k=1)
    m = spcr(x=x[51:250,],y=y1[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    N_1_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[51:250,],y=y1[51:250],k=1,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y1[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    AN_1_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[1:50,], y=y1[1:50], k=10)
    m = spcr(x=x[1:50,],y=y1[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    N_10_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[1:50,],y=y1[1:50],k=10,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y1[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    AN_10_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[51:250,], y=y1[51:250], k=10)
    m = spcr(x=x[51:250,],y=y1[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    N_10_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[51:250,],y=y1[51:250],k=10,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y1[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    AN_10_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[1:50,], y=y2[1:50], k=1)
    m = spcr(x=x[1:50,],y=y2[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    N_1_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[1:50,],y=y2[1:50],k=1,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y2[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    AN_1_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[51:250,], y=y2[51:250], k=1)
    m = spcr(x=x[51:250,],y=y2[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    N_1_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[51:250,],y=y2[51:250],k=1,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y2[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    AN_1_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[1:50,], y=y2[1:50], k=10)
    m = spcr(x=x[1:50,],y=y2[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    N_10_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[1:50,],y=y2[1:50],k=10,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y2[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    AN_10_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[51:250,], y=y2[51:250], k=10)
    m = spcr(x=x[51:250,],y=y2[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    N_10_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
    cv = cv.spcr(x=x[51:250,],y=y2[51:250],k=10,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y2[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/2
    AN_10_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/8
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
ans = data1b()
haha = cbind(mean = ans$mean,sd = ans$sd)
write.table(haha,"data1b.csv")