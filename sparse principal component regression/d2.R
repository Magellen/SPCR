####simulation####
library(mvtnorm)#use this package to generate normal dist
library(spcr)#use the packages to speed up this boring task
options(digits=3)#for simplicity
set.seed(666)
setwd("D:/study/DS2-1/stat760/final")
#case 2
data2 = function(){
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
    nonzero = c(T,F,T,T,F,T,T,F,T,rep(FALSE,11))
    error1 = rnorm(n=1250,mean = 0,sd=0.1)
    error2 = rnorm(n=1250,mean = 0,sd=1)
    sigma21=diag(rep(0,9))
    for(i in 1:9)
    {
      for(j in 1:9)
      {
        sigma21[i,j]=0.9^abs(i-j)
      }
    }
    sigma22=matrix(0,nrow=9,ncol=11)
    sigma23=matrix(0,nrow=11,ncol=9)
    sigma24=diag(rep(1,11))
    sigma2=rbind(cbind(sigma21,sigma22),cbind(sigma23,sigma24))
    x2 = rmvnorm(n=1250,mean = rep(0,20),sigma = sigma2,method="chol")
    beta2=c(-1,0,1,1,0,-1,-1,0,1,rep(0,11))
    y21 = 4*x2%*%beta2+error1
    y22 = 4*x2%*%beta2+error2
    x = x2
    y1 = y21
    y2 = y22
    cv = cv.spcr(x=x[1:50,], y=y1[1:50], k=1)
    m = spcr(x=x[1:50,],y=y1[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    N_1_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[1:50,],y=y1[1:50],k=1,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y1[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    AN_1_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[51:250,], y=y1[51:250], k=1)
    m = spcr(x=x[51:250,],y=y1[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    N_1_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[51:250,],y=y1[51:250],k=1,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y1[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    AN_1_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[1:50,], y=y1[1:50], k=10)
    m = spcr(x=x[1:50,],y=y1[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    N_10_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[1:50,],y=y1[1:50],k=10,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y1[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_50_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_50_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    AN_10_50_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[51:250,], y=y1[51:250], k=10)
    m = spcr(x=x[51:250,],y=y1[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    N_10_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[51:250,],y=y1[51:250],k=10,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y1[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_200_0.1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_200_0.1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    AN_10_200_0.1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[1:50,], y=y2[1:50], k=1)
    m = spcr(x=x[1:50,],y=y2[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    N_1_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[1:50,],y=y2[1:50],k=1,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y2[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    AN_1_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[51:250,], y=y2[51:250], k=1)
    m = spcr(x=x[51:250,],y=y2[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_1_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_1_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    N_1_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[51:250,],y=y2[51:250],k=1,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y2[51:250],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_1_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_1_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    AN_1_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[1:50,], y=y2[1:50], k=10)
    m = spcr(x=x[1:50,],y=y2[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    N_10_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[1:50,],y=y2[1:50],k=10,adaptive=TRUE)
    m = spcr(x=x[1:50,],y=y2[1:50],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_50_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_50_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    AN_10_50_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[51:250,], y=y2[51:250], k=10)
    m = spcr(x=x[51:250,],y=y2[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
    MSE_10_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    P_10_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    N_10_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
    cv = cv.spcr(x=x[51:250,],y=y2[51:250],k=10,adaptive=TRUE)
    m = spcr(x=x[51:250,],y=y2[51:250],k=10,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
    MSEA_10_200_1[t] = mean((y1[251:1250]-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(x[251:1250,]))^2)
    AP_10_200_1[t] = sum(nonzero&(c(m$loadings.B%*%m$gamma)!=0))/6
    AN_10_200_1[t] = sum(!nonzero&(c(m$loadings.B%*%m$gamma)==0))/14
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
ans = data2()
haha = cbind(mean = ans$mean,sd = ans$sd)
write.table(haha,"data2.csv")