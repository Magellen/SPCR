#################################################
####P12 sparse principle component regression####
#################################################

####function####
#1st is soft-thresholding operator
#2nd is original spcr algorithm, the author wrote a R package but in fact he use C++ to achieve it.
#you can see the detail in web: https://cran.r-project.org/web/packages/spcr/index.html
#3rd is adaptive spcr,it first get the parameters in original algorithm and then use these to get a more sparse Beta estimator.
#4th is the function to generate initial lambda, and the method comes from the glmnet packages.
#5th is the function for corss-validation
softsh = function(z,eta){
  sss = 0
  if(eta < abs(z)){
    if(z>0)
    {
      sss = z-eta
      return(sss)
    }else
    {
      sss = z+eta
      return(sss)
    }
  }else
  {
    sss = 0
    return(sss)
  }
}
spcr_yilun = function(x,y,k,lambda_beta,lambda_gamma,xi=0.01,w=0.1,center=TRUE,scale=FALSE){
  if( !is.matrix(x) ) stop("x must be a matrix.")
  if( mode(x)!="numeric" ) stop("x must be numeric.")
  if ( !is.vector(y) ) stop("y must be a vector.")
  if( mode(y)!="numeric" ) stop("y must be numeric.")
  if( center==TRUE ) x <- sweep(x, 2, apply(x,2,mean))
  if( scale==TRUE ) x <- scale(x)
  A <- as.matrix(eigen(var(x))$vectors[ ,1:k])
  gamma0 <- mean(y)
  gamma <- rep(0, k)
  Beta <- matrix( 0, nrow(A), k )
  n = nrow(x)#number of sample
  p = ncol(x)#number of variable 
  k = ncol(A)#number of principle component
  pk = p*k
  pkkone = pk+k+1#number of varialbe in spcr
  #put all the variable in para_new
  para_new = numeric(pkkone)
  for(i in 1:pk){
    para_new[i] = c(Beta)[i]
  }
  for(i in (pk+1):(pk+k)){
    para_new[i] = gamma[i-pk]
  }
  para_new[pk+k+1] = gamma0
  para_diff = rep(10,pkkone)
  #initial value 
  diff_max = 1
  n_itr = 0
  #sum of x^2 dim p*p
  X_inn = t(x)%*%x
  #sum of xy length = p
  XinY = t(x)%*%y
  ####update####
  while(diff_max>0.001){
    para_old = para_new
    ####estimate Beta####
    #y_star dim n*k
    y_star = x%*%A
    #unknown length = n
    unit_vec = rep(gamma0,n)
    #gamma0 * sum of x length = p
    Xingamma0 = t(x)%*%unit_vec
    #sum of x*y_star dim p*k
    XinY_star = t(x)%*%y_star
    #sum of gamma*beta*x^2 length = p 
    X_innBetainGamma = X_inn%*%(Beta%*%gamma)
    #sum of beta*x^2 dim p*k 
    X_innBeta = X_inn%*%Beta
    for(j in 1:k){
      for(l in 1:p){
        residual_A = XinY[l]-Xingamma0[l]-X_innBetainGamma[l]
        residual_B = XinY_star[l,j]-X_innBeta[l,j]
        s = (1-w)*gamma[j]*residual_A+w*residual_B+Beta[l,j]*X_inn[l,l]*(w+(1-w)*gamma[j]^2)
        Beta_old = Beta[l,j]
        Beta[l,j] = softsh(s,0.5*lambda_beta*(1-xi))/(((1-w)*gamma[j]^2+w)*X_inn[l,l]+lambda_beta*xi)
        if(Beta[l,j]!=0){
          for(i in 1:p){
            X_innBetainGamma[i] = X_innBetainGamma[i]+gamma[j]*Beta[l,j]*X_inn[i,l]-gamma[j]*Beta_old*X_inn[i,l]
            X_innBeta[i,j] = X_innBeta[i,j]+Beta[l,j]*X_inn[i,l]-Beta_old*X_inn[i,l]
          }
        }
        if(Beta_old!=0&Beta[l,j]==0){
          for(i in 1:p){
            X_innBetainGamma[i] = X_innBetainGamma[i]-gamma[j]*Beta_old*X_inn[i,l]
            X_innBeta[i,j] = X_innBeta[i,j]-Beta_old*X_inn[i,l]
          }
        }
      }
    }
    ####estimate gamma####
    #x_star dim n*k
    x_star = x%*%Beta
    #sum of square of x_star dim k*k
    X_star_inn = t(x_star)%*%x_star
    #sum of x_star * y length k
    X_starinY = t(x_star)%*%y
    #sum of x_star * gamma0 length k
    X_staringamma0 = t(x_star)%*%unit_vec
    #sum of x_star^2 * gamma length k 
    X_starinngamma = X_star_inn%*%gamma
    for(l in 1:k){
      residual_A = X_starinY[l]-X_staringamma0[l]-X_starinngamma[l]
      s = residual_A+gamma[l]*X_star_inn[l,l]
      gamma_old = gamma[l]
      gamma[l] = softsh((1-w)*s,0.5*lambda_gamma)/((1-w)*X_star_inn[l,l])
      if(is.na(gamma[l])){gamma[l]=0}
      if(gamma[l]!=0){
        for(i in 1:k){
          X_starinngamma[i]=X_starinngamma[i]+gamma[l]*X_star_inn[i,l]-gamma_old*X_star_inn[i,l]
        }
      }
      if(gamma_old!=0&gamma[l]==0){
        for(i in 1:k){
          X_starinngamma[i]=X_starinngamma[i]-gamma_old*X_star_inn[i,l]
        }
      }
    }
    ####estimate gamma0####
    gamma0 = mean(y-x%*%Beta%*%gamma)
    ####estimate A####
    M = svd(t(x)%*%x%*%Beta)
    A = M$u%*%t(M$v)
    #update para_new
    for(i in 1:pk){
      para_new[i] = c(Beta)[i]
    }
    for(i in (pk+1):(pk+k)){
      para_new[i] = gamma[i-pk]
    }
    para_new[pk+k+1] = gamma0
    #update para_diff
    for(i in 1:pkkone){
      para_diff[i] = abs(para_new[i]-para_old[i])
    }
    #update diff_max
    diff_max = max(para_diff)
    if(mean(para_diff)==0)break
    n_itr = n_itr+1
  }
  ans = list(Beta=Beta,gamma=gamma,gamma0=gamma0,A=A,n_itr=n_itr)
  return(ans)
}
adaspcr_yilun = function(x,y,k,lambda_beta,lambda_gamma,xi=0.01,w=0.1,center=TRUE,scale=FALSE){
  spcr.object = spcr_yilun(x,y,k,lambda_beta,lambda_gamma,xi=0.01,w=0.1,center=TRUE,scale=FALSE)
  Beta = spcr.object[[1]]
  gamma = spcr.object[[2]]
  gamma0 = spcr.object[[3]]
  A = spcr.object[[4]]
  if(sum(abs(Beta))==0){
    ans = list(Beta=Beta,gamma=gamma,gamma0=gamma0,A=A,n_itr=0)
    return(ans)
  }
  BetaWeight = Beta/sum(abs(Beta))
  if( center==TRUE ) x <- sweep(x, 2, apply(x,2,mean))
  if( scale==TRUE ) x <- scale(x)
  n = nrow(x)#number of sample
  p = ncol(x)#number of variable 
  k = ncol(A)#number of principle component
  pk = p*k
  pkkone = pk+k+1#number of varialbe in spcr
  #put all the variable in para_new
  para_new = numeric(pkkone)
  for(i in 1:pk){
    para_new[i] = c(Beta)[i]
  }
  for(i in (pk+1):(pk+k)){
    para_new[i] = gamma[i-pk]
  }
  para_new[pk+k+1] = gamma0
  para_diff = rep(10,pkkone)
  #initial value 
  diff_max = 1
  n_itr = 0
  #sum of x^2 dim p*p
  X_inn = t(x)%*%x
  #sum of xy length = p
  XinY = t(x)%*%y
  ####update####
  while(diff_max>0.001){
    para_old = para_new
    ####estimate Beta####
    #y_star dim n*k
    y_star = x%*%A
    #unknown length = n
    unit_vec = rep(gamma0,n)
    #gamma0 * sum of x length = p
    Xingamma0 = t(x)%*%unit_vec
    #sum of x*y_star dim p*k
    XinY_star = t(x)%*%y_star
    #sum of gamma*beta*x^2 length = p 
    X_innBetainGamma = X_inn%*%(Beta%*%gamma)
    #sum of beta*x^2 dim p*k 
    X_innBeta = X_inn%*%Beta
    for(j in 1:k){
      for(l in 1:p){
        residual_A = XinY[l]-Xingamma0[l]-X_innBetainGamma[l]
        residual_B = XinY_star[l,j]-X_innBeta[l,j]
        s = (1-w)*gamma[j]*residual_A+w*residual_B+Beta[l,j]*X_inn[l,l]*(w+(1-w)*gamma[j]^2)
        Beta_old = Beta[l,j]
        Beta[l,j] = softsh(s,0.5*lambda_beta*(1-xi)/(abs(BetaWeight[l,j])+0.0000001))/(((1-w)*gamma[j]^2+w)*X_inn[l,l]+lambda_beta*xi)
        if(Beta[l,j]!=0){
          for(i in 1:p){
            X_innBetainGamma[i] = X_innBetainGamma[i]+gamma[j]*Beta[l,j]*X_inn[i,l]-gamma[j]*Beta_old*X_inn[i,l]
            X_innBeta[i,j] = X_innBeta[i,j]+Beta[l,j]*X_inn[i,l]-Beta_old*X_inn[i,l]
          }
        }
        if(Beta_old!=0&Beta[l,j]==0){
          for(i in 1:p){
            X_innBetainGamma[i] = X_innBetainGamma[i]-gamma[j]*Beta_old*X_inn[i,l]
            X_innBeta[i,j] = X_innBeta[i,j]-Beta_old*X_inn[i,l]
          }
        }
      }
    }
    ####estimate gamma####
    #x_star dim n*k
    x_star = x%*%Beta
    #sum of square of x_star dim k*k
    X_star_inn = t(x_star)%*%x_star
    #sum of x_star * y length k
    X_starinY = t(x_star)%*%y
    #sum of x_star * gamma0 length k
    X_staringamma0 = t(x_star)%*%unit_vec
    #sum of x_star^2 * gamma length k 
    X_starinngamma = X_star_inn%*%gamma
    for(l in 1:k){
      residual_A = X_starinY[l]-X_staringamma0[l]-X_starinngamma[l]
      s = residual_A+gamma[l]*X_star_inn[l,l]
      gamma_old = gamma[l]
      gamma[l] = softsh((1-w)*s,0.5*lambda_gamma)/((1-w)*X_star_inn[l,l])
      if(is.na(gamma[l])){gamma[l]=0}
      if(gamma[l]!=0){
        for(i in 1:k){
          X_starinngamma[i]=X_starinngamma[i]+gamma[l]*X_star_inn[i,l]-gamma_old*X_star_inn[i,l]
        }
      }
      if(gamma_old!=0&gamma[l]==0){
        for(i in 1:k){
          X_starinngamma[i]=X_starinngamma[i]-gamma_old*X_star_inn[i,l]
        }
      }
    }
    ####estimate gamma0####
    gamma0 = mean(y-x%*%Beta%*%gamma)
    ####estimate A####
    M = svd(t(x)%*%x%*%Beta)
    A = M$u%*%t(M$v)
    #update para_new
    for(i in 1:pk){
      para_new[i] = c(Beta)[i]
    }
    for(i in (pk+1):(pk+k)){
      para_new[i] = gamma[i-pk]
    }
    para_new[pk+k+1] = gamma0
    #update para_diff
    for(i in 1:pkkone){
      para_diff[i] = abs(para_new[i]-para_old[i])
    }
    #update diff_max
    diff_max = max(para_diff)
    if(mean(para_diff)==0)break
    n_itr = n_itr+1
  }
  ans = list(Beta=Beta,gamma=gamma,gamma0=gamma0,A=A,n_itr=n_itr)
  return(ans)
}
ini.lambda = function(x, y, k, w, xi){
  A <- as.matrix(eigen(var(x))$vectors[ ,1:k])
  y_star <- x %*% A
  lambda_beta_max <- 0
  for( j in 1:ncol(y_star) )
  {
    for( l in 1:ncol(x) )
    {
      s <- t(x[ , l ]) %*% y_star[ , j ] * 2 * w / ( 1 - xi )
      if( lambda_beta_max < abs( s ) ) lambda_beta_max <- abs( s )
    }
  }
  ss <- lambda_beta_max - lambda_beta_max*10^(-5)
  return(ss)
}
cv_spcr_yilun <- function(x, y, k, w=0.1, xi=0.01, nfolds=5, adaptive=FALSE, center=TRUE, scale=FALSE, 
                          lambda_beta_length=10, lambda_gamma_length=10, lambda_beta=NULL, lambda_gamma=NULL){
  if( !is.matrix(x) ) stop("x must be a matrix.")
  if( mode(x)!="numeric" ) stop("x must be numeric.")
  if ( !is.vector(y) ) stop("y must be a vector.")
  if( mode(y)!="numeric" ) stop("y must be numeric.")
  if( k < 1 ) stop("k is an integer and larger than one.")
  n = length(y)
  #get initial lambda by glm algorithm
  ini.lambda.beta = ini.lambda.gamma = ini.lambda( x=x, y=y, k=k, w=w, xi=xi )
  #get a path of lambda
  lambda.beta.candidate = rev( seq( n*0.005, c(ini.lambda.beta), length=lambda_beta_length ) )
  lambda.gamma.candidate = rev( seq(n* 0.005, c(ini.lambda.gamma), length=lambda_gamma_length ) )
  #use path if provided
  if( is.null(lambda_beta) != TRUE ) lambda.beta.candidate = sort(lambda_beta, decreasing=TRUE)
  if( is.null(lambda_gamma) != TRUE ) lambda.gamma.candidate = sort(lambda_gamma, decreasing=TRUE)
  #initial parameters
  A.ini <- as.matrix(eigen(var(x))$vectors[ ,1:k])
  gamma0.ini <- mean(y)
  gamma.ini <- rep(0, k)
  Beta.ini <- matrix( 0, nrow(A.ini), k )
  CV.mat = matrix( 0, length(lambda.gamma.candidate), length(lambda.beta.candidate) )
  foldid = sample(rep(seq(nfolds),length=n))
  x.all = x
  y.all = y
  ####start cross validation
  for(i in seq(nfolds))
  {
    num.foldid <- which(foldid==i)
    x <- x.all[ -num.foldid, ]
    y <- y.all[ -num.foldid ]
    x.test.cv <- x.all[ num.foldid, ]
    y.test.cv <- y.all[ num.foldid ]
    #center or scale
    if( center==TRUE ){
      x_ori  <- x
      x <- sweep(x_ori, 2, apply(x_ori,2,mean))
      x.test.cv <- sweep(x.test.cv, 2, apply(x_ori,2,mean))
    }
    if( scale==TRUE ){
      x_ori  <- x
      x <- scale(x_ori)
      x.test.cv <- sweep(sweep(x.test.cv, 2, apply(x_ori, 2, mean)), 2, apply(x_ori, 2, sd), FUN="/")
    }
    for( itr.lambda.gamma in 1:length(lambda.gamma.candidate) )
    {
      lambda.gamma <- lambda.gamma.candidate[itr.lambda.gamma]
      A <- A.ini
      gamma0 <- gamma0.ini
      gamma <- gamma.ini
      Beta <- Beta.ini
      for( itr.lambda.beta in 1:length(lambda.beta.candidate) )
      {
        lambda.beta <- lambda.beta.candidate[itr.lambda.beta]
        para_old <- c(gamma0, gamma, matrix(Beta, 1, nrow(Beta)*ncol(Beta)))
        para_new <- para_old + 10
        if( adaptive==FALSE ){
          spcr.object <- spcr_yilun(x,y,k,lambda_beta=lambda.beta,lambda_gamma=lambda.gamma,xi,w)
          Beta <- spcr.object[[1]]
          gamma <- spcr.object[[2]]
          gamma0 <- spcr.object[[3]]
          A <- spcr.object[[4]]
        } else {
          spcr.object <- adaspcr_yilun(x,y,k,lambda_beta=lambda.beta,lambda_gamma=lambda.gamma,xi,w)
          Beta <- spcr.object[[1]]
          gamma <- spcr.object[[2]]
          gamma0 <- spcr.object[[3]]
          A <- spcr.object[[4]]
        }
        para_old <- para_new
        para_new <- c(gamma0, gamma, c(Beta))
        if( mean(abs(para_new-para_old)) == 0 ) break
        #CV-error
        s_cv <- mean( ( y.test.cv - gamma0 - t(gamma) %*% t(Beta) %*% t(x.test.cv) )^2 )
        
        #Strock of CV-error
        CV.mat[ itr.lambda.gamma, itr.lambda.beta ] <- CV.mat[ itr.lambda.gamma, itr.lambda.beta ] + s_cv/nrow(x.test.cv)				
      }
    }
  }
  CV.mat <- CV.mat/nfolds
  ####Search of min CV
  minCandi.col <- whichiminCandi.col <- rep(0, nrow(CV.mat))
  for(i in 1:nrow(CV.mat))
  {
    whichiminCandi.col[i] <- which.min(CV.mat[i, ])
    minCandi.col[i] <- min(CV.mat[i, ])
  }
  whichiminCandi.row <- which.min( CV.mat[ , whichiminCandi.col[ which.min(minCandi.col) ]] )
  minCandi.row <- min( CV.mat[ , whichiminCandi.col[ which.min(minCandi.col) ]] )		
  ####Selected tuning parameters by CV
  lambda.gamma.cv <- lambda.gamma.candidate[ whichiminCandi.row ]
  lambda.beta.cv <- lambda.beta.candidate[ whichiminCandi.col[ which.min(minCandi.col) ] ]
  
  ans <- list( lambda.gamma.seq=lambda.gamma.candidate, lambda.B.seq=lambda.beta.candidate, CV.mat=CV.mat, lambda.gamma.cv=lambda.gamma.cv, lambda.B.cv=lambda.beta.cv, cvm=min(minCandi.row), call=match.call() )
  return(ans)
} 

####test the custom function and R package####
library(mvtnorm)
library(spcr)
set.seed(666)
x1a[1:250,] = rmvnorm(n=250,mean = rep(0,10),sigma = diag(rep(1,10)),method="chol")
beta1a = c(2,1,rep(0,8))
error1 = rnorm(n=1250,mean = 0,sd=0.1)
error2 = rnorm(n=1250,mean = 0,sd=1)
y1a1 = x1a%*%beta1a+error1
y1a2 = x1a%*%beta1a+error2
x = x1a
y1 = y1a1
y2 = y1a2
#spcr package
set.seed(666)
cv = cv.spcr(x=x[1:50,], y=y1[1:50], k=1)
m = spcr(x=x[1:50,],y=y1[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv)
#spcr custom
set.seed(666)
cv1 = cv_spcr_yilun(x=x[1:50,], y=y1[1:50], k=1)
m1 = spcr_yilun(x=x[1:50,],y=y1[1:50],k=1,lambda_beta=cv$lambda.B.cv,lambda_gamma=cv$lambda.gamma.cv)
#adaspcr package
set.seed(666)
cv = cv.spcr(x=x[1:50,], y=y1[1:50], k=1,adaptive = TRUE)
m = spcr(x=x[1:50,],y=y1[1:50],k=1,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE)
#adaspcr custom
set.seed(666)
cv1 = cv_spcr_yilun(x=x[1:50,], y=y1[1:50], k=1,adaptive = TRUE)
m1 = adaspcr_yilun(x=x[1:50,],y=y1[1:50],k=1,lambda_beta=cv$lambda.B.cv,lambda_gamma=cv$lambda.gamma.cv)
#package time for cross validation (default 100 time search for the best hyper parameters.)
set.seed(666)
system.time(cv.spcr(x=x[1:50,], y=y1[1:50], k=1,adaptive = TRUE))#0.31s
#custom time for cross validation (default 100 time search for the best hyper parameters.)
set.seed(666)
system.time(cv_spcr_yilun(x=x[1:50,], y=y1[1:50], k=1,adaptive = TRUE))#10.6s

#the parameters B, gamma, gamma0 are all the same but running C++ is about 30 time faster than R. 
#It will cost you several days to run the simulation part and actually the numbers in those table 
#are highly effected by the random seed. If you use matlab you could follow the path of my custom 
#function because matlab is comparable to C++ and but if you want to use python then you had better to
#import the C++ module to speed up. At the beginning I only use the custom function but my computer
#cannot finish the simulation for 1 table after 2 days calculation, so I switch to R packages for
#simulation. I can give you the code for custom function simulation but you will not want to run it again.
#There are 6 files you need to run. d1a, d1b, d2, d3a, d3b, and real_data. You can use remote server or 
#open 6 Rstudio to run the code at the same time. There are totally 6 tables in the paper. When you finish
#the first 5 files you can get the first 4 table. Table 5 just a summary of the real data so I skip it.
#For table 6, since the housing dataset is no longer exist in the UCI website so I only run the code 
#based on the 4 database. For these 6 files, I use the R packages to speed up since the calculation in R is 
#too slow. I almost follow the algorithm describe in the C++ and there is no much space to optimize the code,
#so I think my function is equivalent to the package. The only difference is just the platform and I am very
#upset about the R language now :)