rm(list = ls())
set.seed(8123)
library(randomForest)
library(SuperLearner)
source('data_prep.R') #X has all 71 variables
source('toSource.R')
data_prep(d=71)
data <- list(A=matrix(A), Y=matrix(Y), Z=Z, W=W, X=X)
N = length(Y)
n = floor(N/2)
index_train = sample(N,n)
index_eval = seq(N)[-index_train]
X = X[index_train,]; Y = Y[index_train]; A = A[index_train]
random_sample = list(index_train, index_eval)
method_pi = c("SL.mean","SL.randomForest","SL.glm")
method_mu = c("SL.mean","SL.randomForest","SL.glm")
pi_sl = function(xnew,ytrain,xtrain,method=method_pi){
  #bw <- npcdistbw(r_train0~x_train0)
  xnew = data.frame(xnew)
  sl_pi = SuperLearner(Y = ytrain, X = data.frame(xtrain), family = binomial(),
                       SL.library = c("SL.mean","SL.glm","SL.bartMachine","SL.ridge"))
  pred = predict(sl_pi, xnew, onlySL = TRUE)
  return(pred$pred[, 1])
}
pihat_sl = pi_sl(data$X[-index_train,],as.numeric(A)-1,X)
summary(pihat_sl)

mu_sl=function(xnew,ytrain,xtrain,method=method_mu){
  xnew = data.frame(xnew)
  sl_m = SuperLearner(Y = ytrain, X = data.frame(xtrain),
                      SL.library = method)
  pred = predict(sl_m, xnew, onlySL = TRUE)
  return(pred$pred[,1])
}

mu1hat_sl = mu_sl(data$X[-index_train,],as.numeric(Y[A==1]),X[A==1,])
mu0hat_sl = mu_sl(data$X[-index_train,],as.numeric(Y[A==0]),X[A==0,])

A = data$A[-index_train]
X = data$X[-index_train,]
Y = data$Y[-index_train]
pseudo_sl <- ((A-pihat_sl)/(pihat_sl*(1-pihat_sl)))*(Y-A*mu1hat_sl-(1-A)*mu0hat_sl) + mu1hat_sl-mu0hat_sl
all = list(pseudo_sl,X,A,Y,random_sample)
save(all, file = "pseudo_unconfound.RData")

