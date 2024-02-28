rm(list = ls())
set.seed(8123)
library(caret)
library(corpcor)
library(rdetools)
library(expm)
library(boot)
#library(SuperLearner)
#library(gam)
#library(polspline)
library(ranger)
library(glmnet)
#library(SIS)
#library(tmle)
#library(twang)
library(gbm)
source('data_prep.R') #X has all 11 variables
source('toSource.R')
CF_K_fold = 2
lm_Hh1 <- 0.001
lm_Qh1 <- 0.01
lm_Hh0 <- 0.001
lm_Qh0 <- 0.01
lm_Hq1 <- 0.01
lm_Qq1 <- 0.001
lm_Hq0 <- 0.01
lm_Qq0 <- 0.001
kernel_sigma_h <- 35
kernel_sigma_q <- 20
data_prep(d=7)
data <- list(A=matrix(A), Y=matrix(Y), Z=Z, W=W, X=X)
random_sample <- createFolds(1:nrow(data$Y), k = CF_K_fold)
data1 <- sapply(data, function(x) matrix(x[random_sample[[1]], ], ncol = ncol(x)))
eval_data <- sapply(data, function(x) matrix(x[random_sample[[2]], ], ncol = ncol(x)))
est_proxy <- CF_cate(data1, kernel_sigma_h, kernel_sigma_q, lm_Hh1, lm_Qh1, lm_Hh0, lm_Qh0, lm_Hq1, lm_Qq1, lm_Hq0, lm_Qq0, CF_K_fold, eval_data = eval_data)
pseudo_proxy = est_proxy[[1]]
all = list(pseudo_proxy,X[random_sample[[2]], ],A[random_sample[[2]]],Y[random_sample[[2]] ],Z[random_sample[[2]], ],W[random_sample[[2]], ],random_sample)
save(all, file = "pseudo_proxy.RData")

