rm(list = ls())
set.seed(1234)
library(MASS)
library(splines)
source('toSource.R')
load("pseudo_proxy.RData")
pseudo_proxy <- all[[1]]
random_sample <- all[[length(all)]]
X <- all[[2]]
A <- all[[3]]
Y <- all[[4]]
X1 = X[,"surv2md1"]


k = 0
D1 = D2 = D3 = 0
x1_list = seq(min(X1)+0.1,max(X1)-0.01,by = 0.1)
l = length(x1_list)
sigma = sd(pseudo_proxy)
est_forster_bs = sd_forster_bs = rep(NA,l)
est_forster_ns = sd_forster_ns = rep(NA,l)
est_forster_poly = sd_forster_poly = rep(NA,l)
x1 = x1_list[1]
for (x1 in x1_list){
  k = k+1
  cat(paste0(k,"/", length(x1_list),' '))
  tau_forster_bs = series_cv_new(X1,pseudo_proxy,x1,  type = "forster", basis_type = "bs", KK=5)
  #g(X, Y, x_pred, df,type, basis_type)  %=%  list(X1, pseudo_proxy, x1, 50, "forster", "bs")
  est_forster_bs[k] <- tau_forster_bs[[1]]
  D1_list = tau_forster_bs[[3]]
  var1_list = sapply(D1_list, function(D_1){ 
    tmp = series_df(X1,pseudo_proxy,x1 ,df=D_1,type = "forster", basis_type = "bs", std=TRUE)
    tmp[[2]]
    })
  D1 = c(D1,D1_list)
  sd_forster_bs[k] = mean(var1_list)
  
  
  tau_forster_poly = series_cv_new(X1,pseudo_proxy,x1,  type = "forster", basis_type = "poly", KK=5)
  #g(X, Y, x_pred, df,type, basis_type)  %=%  list(X1, pseudo_proxy, x1, 5, "forster", "poly")
  est_forster_poly[k] <- tau_forster_poly[[1]]
  D3_list = tau_forster_poly[[3]]
  var3_list = sapply(D3_list, function(D_3){ 
    tmp = series_df(X1,pseudo_proxy,x1 ,df=D_3,type = "forster", basis_type = "poly",std=TRUE)
    tmp[[2]]
  })
  D3 = c(D3,D3_list)
  sd_forster_poly[k] = mean(var3_list)
}

proxy_data = list(est_forster_poly,sd_forster_poly,est_forster_bs,sd_forster_bs,x1_list)
save(proxy_data, file = "proxy_plotdata.RData")



