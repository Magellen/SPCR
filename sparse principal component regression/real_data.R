rm(list=ls())
library(xlsx)
library(spcr)
####real data####
#I dont find the housing data in the UCI database.
setwd("D:/study/DS2-1/stat760/final/real data")#change this line to your directory.

####energy data####
energy = read.xlsx("ENB2012_data.xlsx",1)
energy = data.frame(energy)
energy$NA.=NULL
energy = energy[complete.cases(energy),]
energy[,1:8] = scale(energy[,1:8])
#dim 1296 * 10
SPCR1 = numeric(100)
aSPCR1 = numeric(100)
SPCR2 = numeric(100)
aSPCR2 = numeric(100)
set.seed(666)
for(t in 1:100){
  print(t)
  sam = sample(1:nrow(energy),100)
  train = energy[sam,]
  test = energy[-sam,]
  cv = cv.spcr(x=as.matrix(train[,1:8]), y=train$Y1, k=5,center = FALSE)
  m = spcr(x=as.matrix(train[,1:8]), y=train$Y1, k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,center = FALSE)
  SPCR1[t] = mean((test$Y1-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test[,1:8]))^2)
  cv = cv.spcr(x=as.matrix(train[,1:8]), y=train$Y1, k=5,adaptive = TRUE,center = FALSE)
  m = spcr(x=as.matrix(train[,1:8]), y=train$Y1, k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE,center = FALSE)
  aSPCR1[t] = mean((test$Y1-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test[,1:8]))^2)
  cv = cv.spcr(x=as.matrix(train[,1:8]), y=train$Y2, k=5,center = FALSE)
  m = spcr(x=as.matrix(train[,1:8]), y=train$Y2, k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,center = FALSE)
  SPCR2[t] = mean((test$Y2-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test[,1:8]))^2)
  cv = cv.spcr(x=as.matrix(train[,1:8]), y=train$Y2, k=5,adaptive = TRUE,center = FALSE)
  m = spcr(x=as.matrix(train[,1:8]), y=train$Y2, k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE,center = FALSE)
  aSPCR2[t] = mean((test$Y2-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test[,1:8]))^2)
}

####forestfires####
fire = read.csv("forestfires.csv")
fire$month = NULL
fire$day = NULL
fire[,1:10] = scale(fire[,1:10])
#dim 517 11
SPCR3 = numeric(100)
aSPCR3 = numeric(100)
set.seed(666)
for(t in 1:100){
  print(t)
  sam = sample(1:nrow(fire),100)
  train = fire[sam,]
  test = fire[-sam,]
  cv = cv.spcr(x=as.matrix(train[,1:10]), y=train$area, k=5,center = FALSE)
  m = spcr(x=as.matrix(train[,1:10]), y=train$area, k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,center = FALSE)
  SPCR3[t] = mean((test$area-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test[,1:10]))^2)
  cv = cv.spcr(x=as.matrix(train[,1:10]), y=train$area, k=5,adaptive = TRUE,center = FALSE)
  m = spcr(x=as.matrix(train[,1:10]), y=train$area, k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE,center = FALSE)
  aSPCR3[t] = mean((test$area-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test[,1:10]))^2)
}

####concrete####
con = read.xlsx("Concrete_Data.xls",1)
con[,1:8] = scale(con[,1:8])
#dim 1030 9
SPCR4 = numeric(100)
aSPCR4 = numeric(100)
set.seed(666)
for(t in 1:100){
  print(t)
  sam = sample(1:nrow(con),100)
  train = con[sam,]
  test = con[-sam,]
  cv = cv.spcr(x=as.matrix(train[,1:8]), y=train$Concrete.compressive.strength.MPa..megapascals.., k=5,center = FALSE)
  m = spcr(x=as.matrix(train[,1:8]), y=train$Concrete.compressive.strength.MPa..megapascals.., k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,center = FALSE)
  SPCR4[t] = mean((test$Concrete.compressive.strength.MPa..megapascals..-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test[,1:8]))^2)
  cv = cv.spcr(x=as.matrix(train[,1:8]), y=train$Concrete.compressive.strength.MPa..megapascals.., k=5,adaptive = TRUE,center = FALSE)
  m = spcr(x=as.matrix(train[,1:8]), y=train$Concrete.compressive.strength.MPa..megapascals.., k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE,center = FALSE)
  aSPCR4[t] = mean((test$Concrete.compressive.strength.MPa..megapascals..-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test[,1:8]))^2)
}

####communities####
crime = read.csv("communities.data", header = FALSE, sep = ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "",
                 na.strings="?",strip.white=TRUE,stringsAsFactors = default.stringsAsFactors())
crime = crime[complete.cases(crime),]
crime$V4 = NULL
crime[,1:126] = scale(crime[,1:126])
#dim 123 127
SPCR5 = numeric(100)
aSPCR5 = numeric(100)
SPCR6 = numeric(100)
aSPCR6 = numeric(100)
set.seed(666)
for(t in 1:100){
  print(t)
  sam0 = sample(1:nrow(crime),50)
  train0 = crime[sam0,]
  test0 = crime[-sam0,]
  sam = sample(1:nrow(crime),100)
  train = crime[sam,]
  test = crime[-sam,]
  cv = cv.spcr(x=as.matrix(train0[,1:126]), y=train0$V128, k=5,center = FALSE)
  m = spcr(x=as.matrix(train0[,1:126]), y=train0$V128, k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,center = FALSE)
  SPCR5[t] = mean((test0$V128-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test0[,1:126]))^2)
  cv = cv.spcr(x=as.matrix(train0[,1:126]), y=train0$V128, k=5,adaptive = TRUE,center = FALSE)
  m = spcr(x=as.matrix(train0[,1:126]), y=train0$V128, k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE,center = FALSE)
  aSPCR5[t] = mean((test0$V128-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test0[,1:126]))^2)
  cv = cv.spcr(x=as.matrix(train[,1:126]), y=train$V128, k=5,center = FALSE)
  m = spcr(x=as.matrix(train[,1:126]), y=train$V128, k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,center = FALSE)
  SPCR6[t] = mean((test$V128-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test[,1:126]))^2)
  cv = cv.spcr(x=as.matrix(train[,1:126]), y=train$V128, k=5,adaptive = TRUE,center = FALSE)
  m = spcr(x=as.matrix(train[,1:126]), y=train$V128, k=5,lambda.B=cv$lambda.B.cv,lambda.gamma=cv$lambda.gamma.cv,adaptive = TRUE,center = FALSE)
  aSPCR6[t] = mean((test$V128-m$gamma0-t(m$gamma)%*%t(m$loadings.B)%*%t(test[,1:126]))^2)
}

result = rbind(SPCR1,aSPCR1,SPCR2,aSPCR2,SPCR3,aSPCR3,SPCR4,aSPCR4,SPCR5,aSPCR5,SPCR6,aSPCR6)
MEAN = apply(result,1,mean)
SD = apply(result,1,sd)
ans = list(mean = MEAN,sd = SD)
haha = cbind(mean = ans$mean,sd = ans$sd)
write.table(haha,"realdata.csv")