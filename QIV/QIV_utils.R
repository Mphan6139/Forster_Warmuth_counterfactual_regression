### Original Code From Jiewen

library(transport)
library(ggplot2)
library(MASS)
library(latex2exp)
library(SuperLearner)
library(splines)
library(MASS)
library(mgcv)
source("cust_MySL.R")
source("0.Functions_PMMR.R")
source("cust_fun.R")

construct_pseudo_inv_delta_a = function(a,z,p1,p0,pi){
  pz = p1*z + p0*(1-z)
  return(
    ((2*z-1)/(pi*z+(1-pi)*(1-z))*(a - pz) + p1 - p0) 
    / ((p1 - p0)^2)
  )
}

construct_pseudo_wald = function(y,a,z,p1,p0,e1,e0){
  pz = p1*z + p0*(1-z)
  ez = e1*z + e0*(1-z)
  rho = p1*pi + p0*(1-pi)
  delta_a_inv = 1 / (p1-p0)
  delta = (e1-e0)/(p1-p0)
  
  res = (2*z-1)/(pi*z+(1-pi)*(1-z)) *
    delta_a_inv *
    ( ((y*(1-a)) - ez) -
        (a - pz)*delta ) +
    delta
  return(res)
}

# Y added 
IF = function(y,a,z,p1,p0,e1,e0,pi){
  pz = p1*z + p0*(1-z)
  ez = e1*z + e0*(1-z)
  rho = p1*pi + p0*(1-pi)
  delta_a_inv = 1 / (p1-p0)
  delta = (e1-e0)/(p1-p0)
  
  res = rho *
    (2*z-1)/(pi*z+(1-pi)*(1-z)) *
    delta_a_inv *
    ( ((y*(1-a)) - ez) - 
        (a - pz)*delta ) +
    a*(y + delta)
  return(res)
}

IF_fw = function(y,a,z,p1,p0,e1,e0,delta_a_inv,delta_v){
  pz = p1*z + p0*(1-z)
  ez = e1*z + e0*(1-z)
  rho = p1*pi + p0*(1-pi)
  # delta_v = (e1-e0)/(p1-p0)
  
  res = rho *
    (2*z-1)/(pi*z+(1-pi)*(1-z)) *
    delta_a_inv *
    ( ((y*(1-a)) - ez) -
        (a - pz)*delta_v) +
    a*(y+delta_v)
  return(res)
}

IF_fw2 = function(y,a,z,p1,p0,e1,e0,delta_a_inv,delta_v){
  pz = p1*z + p0*(1-z)
  ez = e1*z + e0*(1-z)
  rho = p1*pi + p0*(1-pi)
  delta_v = (e1-e0)/(p1-p0)
  
  res = rho *
    (2*z-1)/(pi*z+(1-pi)*(1-z)) *
    delta_a_inv *
    ( ((y*(1-a)) - ez) -
        (a - pz)*delta_v) +
    a*(y+delta_v)
  return(res)
}

# super_learner_packed
my_SL = function(data_train, data_pred, choice = c(1,7)){
  SL.hpara <- list()
  SL.hpara$SLL <- choice
  # Superlearner basic learning algorithms:
  # 1: GLM
  # 2: lasso/ridge
  # 3: earth
  # 4: GAM
  # 5: xgboost
  # 6: polynomial spline
  # 7: random forest
  # 9: gbm
  # 10: 1-layer MLP
  SL.hpara$MLPL <- c(1)
  SL.hpara$MTRY <- c(1)
  SL.hpara$NMN <- 50
  SL.hpara$MLPdecay <- 10^c(-1,-3)
  
  ## Estimate pi P(Z=1|X)
  pi.Fit <- MySL(Data = data_train,                # training dataset
                 locY = 4,                         # response variable Z = column 4
                 locX = c(1,2),                    # explanatory variable X1,X2
                 Ydist = binomial(),               # Z is binary
                 SL.list = SL.hpara$SLL,           # Machine learning algorithms
                 MTRY = SL.hpara$MTRY,        
                 MLPL = SL.hpara$MLPL,
                 NMN = SL.hpara$NMN,
                 MLPdecay = SL.hpara$MLPdecay)
  
  ## Estimate p_z P(A=1|Z=z,X)
  pz.Fit <- MySL(Data = data_train,                # training dataset
                 locY = 3,                         # response variable A = column 3
                 locX = c(1,2,4),                  # explanatory variable X1,X2,Z
                 Ydist = binomial(),               # A is binary
                 SL.list = SL.hpara$SLL,           # Machine learning algorithms
                 MTRY = SL.hpara$MTRY,        
                 MLPL = SL.hpara$MLPL,
                 NMN = SL.hpara$NMN,
                 MLPdecay = SL.hpara$MLPdecay)
  
  ## Estimate e_z E[Y(1-A)|Z=z,X]
  ez.Fit <- MySL(Data = data_train,                # training dataset
                 locY = 6,                         # response variable Y = column 1
                 locX = c(1,2,4),                  # explanatory variable X1,X2,Z
                 Ydist = gaussian(),               # Y is continuous (use gaussian in this case)
                 SL.list = SL.hpara$SLL,           # Machine learning algorithms
                 MTRY = SL.hpara$MTRY,        
                 MLPL = SL.hpara$MLPL,
                 NMN = SL.hpara$NMN,
                 MLPdecay = SL.hpara$MLPdecay)
  
  pi = predict(pi.Fit,newdata=data_pred[,1:2])$pred
  p1 = predict(pz.Fit,newdata=data.frame(data_pred[,1:2],z=1))$pred
  p0 = predict(pz.Fit,newdata=data.frame(data_pred[,1:2],z=0))$pred
  e1 = predict(ez.Fit,newdata=data.frame(data_pred[,1:2],z=1))$pred
  e0 = predict(ez.Fit,newdata=data.frame(data_pred[,1:2],z=0))$pred
  
  nuisance_df = data.frame(pi,p1,p0,e1,e0)
  
  return(nuisance_df)
}

my_fw = function(data_train, data_pred, bs = "bs", k = 4, choice = c(1,7), CV = F, df_grid = NULL){
  n_train = dim(data_train)[1]
  n_pred = dim(data_pred)[1]
  s = rep(1:2, each=n_train/2)
  
  df_nuisance_1 = my_fw_aux(data_train = data_train[s==1,], aux = data_train[s==2,],
                            data_pred = data_pred, CV = CV, k = k, bs=bs, choice = choice)
  df_nuisance_2 = my_fw_aux(data_train = data_train[s==2,], aux = data_train[s==1,],
                            data_pred = data_pred, CV = CV, k = k, bs=bs, choice = choice)
  df_nuisance = (df_nuisance_1 + df_nuisance_2)/2
  return(df_nuisance)
}


my_fw_aux = function(data_train, aux, data_pred, bs = "bs", k = 4, choice = c(1,7), CV = F, df_grid = NULL){
  # # TEST
  # data_train = Data[s==1,]
  # aux = Data[s==2,]
  # data_pred = Data[s==3,]
  # choice = c(1,7)
  # CV=T
  # CV=F
  # bs = 'bs'
  # df_grid = NULL
  
  n_train = dim(data_train)[1]
  n_aux = dim(aux)[1]
  n_pred = dim(data_pred)[1]
  
  # use nuisance_df1 to construct the pseudo outcome
  SL_res = my_SL(data_train, rbind(aux,data_pred), choice = choice)
  df_aux = SL_res[1:n_aux,]
  df_pred = SL_res[ (n_aux+1):(n_aux+n_pred), ]
  
  pseudo_pdiff_inv_if = construct_pseudo_inv_delta_a(aux$a,
                                                     aux$z,
                                                     df_aux$p1,
                                                     df_aux$p0,
                                                     df_aux$pi)
  
  pseudo_wald_if = construct_pseudo_wald(aux$y,aux$a,aux$z,
                                         df_aux$p1,df_aux$p0,
                                         df_aux$e1,df_aux$e0)
  
  if(CV == T){
    # Get the optimal basis order from cross fitting for pseudo_pdiff_inv_if1
    cv_res = series_cv_cust(data_pred, pseudo_pdiff_inv_if, type = "forster", basis_type = bs, df_grid = df_grid)
    k_from_cv = cv_res[[1]][1,1:3]
    basis_train = create_basis1(aux$x1,aux$x2, k=k_from_cv, basis_type = bs)
    basis_pred = create_basis1(data_pred$x1,data_pred$x2, k=k_from_cv, basis_type = bs)
  }else{
    basis_train = create_basis1(aux$x1,aux$x2, k=k, basis_type = bs)
    basis_pred = create_basis1(data_pred$x1,data_pred$x2, k=k, basis_type = bs)
  }
  delta_a_inv_pred = series_df_cust(basis_train, pseudo_pdiff_inv_if, basis_pred, type = "forster")[[1]]
  
  if(CV == T){
    # Get the optimal basis order from cross fitting for pseudo_wald_if1
    cv_res = series_cv_cust(data_pred, pseudo_wald_if, type = "forster", basis_type = bs, df_grid = df_grid)
    k_from_cv = cv_res[[1]][1,1:3]
    basis_train = create_basis1(aux$x1,aux$x2, k=k_from_cv, basis_type = bs)
    basis_pred = create_basis1(data_pred$x1,data_pred$x2, k=k_from_cv, basis_type = bs)
  }else{
    basis_train = create_basis1(aux$x1,aux$x2, k=k, basis_type = bs)
    basis_pred = create_basis1(data_pred$x1,data_pred$x2, k=k, basis_type = bs)
  }
  delta_pred = series_df_cust(basis_train, pseudo_wald_if, basis_pred, type = "forster")[[1]]
  
  pi = df_pred$pi
  p1 = df_pred$p1
  p0 = df_pred$p0
  e1 = df_pred$e1
  e0 = df_pred$e0
  
  nuisance_df = data.frame(pi,p1,p0,e1,e0,delta_a_inv_pred,delta_pred)
  return(nuisance_df)
}

two_fold_eval = function(data_input, CV = F, k=4, bs = 'bs', choice = c(1,7)){
  y = data_input$y
  a = data_input$a
  z = data_input$z
  x1 = data_input$x1
  x2 = data_input$x2
  s = data_input$s
  
  Data = data.frame(x1,x2,a,z,y,y*(1-a),s)
  
  get_one = function(combination){
    
    # Test
    # combination = c(1,2)
    # choice = c(1,7)
    # CV = T
    # bs = 'bs'
    # k = 4
    
    train = combination[1]
    pred = combination[2]
    
    nuisance_df_one_step = my_SL(data_train = Data[s==train,],
                                 data_pred = Data[s==pred,], choice = choice)
    nuisance_df_two_step = my_fw(data_train = Data[s==train,],
                                 data_pred = Data[s==pred,], CV = CV, k = k, bs=bs, choice = choice)
    
    estimate_wald = with(nuisance_df_one_step, (y[s==pred] +  (e1-e0)/(p1-p0))
                         * a[s==pred] )
    
    estimate_IF = with(nuisance_df_one_step, IF(y[s==pred],a[s==pred],z[s==pred],
                                                p1,p0,e1,e0))
    estimate_IF_fw = with(nuisance_df_two_step, IF_fw(y[s==pred],a[s==pred],z[s==pred],
                                                      p1,p0,e1,e0,delta_a_inv_pred, delta_pred))
    # estimate_IF_fw = 0
    return(data.frame(estimate_wald,
                      estimate_wald_fw,
                      estimate_IF,
                      estimate_IF_fw))
    
  }
  pred_1 = get_one(c(1,2))
  pred_2 = get_one(c(2,1))
  
  pred_res = rbind(pred_2,pred_1)
  return(pred_res)
}


var_compute = function(estimate, a, s, nc){
  psi_estimate = mean(estimate)
  each_fold_var = rep(NA,nc)
  for(i in 1:nc){
    idx = which(s==i)
    cur_a = a[idx]
    cur_estimate = estimate[idx]
    each_fold_var[i] = sum( ((cur_estimate - cur_a*psi_estimate))^2 ) / length(idx)
  }
  return(mean(each_fold_var))
}

var_compute_m = function(bias_list,var_list){
  tmp = do.call(rbind, bias_list)
  median_theta = apply(tmp, 2, median, na.rm = TRUE)
  adj_var = mapply(FUN=function(X,Y){(X-Y)^2}, X=tmp, Y=median_theta)
  median_var = do.call(rbind, var_list) + adj_var
  median_var = apply(median_var, 2, median, na.rm = TRUE)
  return(median_var)
}

# -----------------------------------------------------------------------------

'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

library(MASS)
library(mgcv)

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

create_basis1 = function(x1,x2,k = 3, basis_type = "bs"){
  
  if (length(k)>1) {
    c1 = as.numeric(k[1])
    c2 = as.numeric(k[2])
    c3 = as.numeric(k[3])
  } else if (length(k) == 1) {
    c1 = c2 = c3 = k
  } else {
    print("The parameter is neither a vector nor a single number.")
  }
  
  if(basis_type == "bs") {
    b1 = bs(x1, df = c1)
    b2 = bs(x2, df = c2)
    b3 = bs(x1*x2, df = c3)
    res = cbind(1, b1, b2, b3)
  }
  else{
    b1 = poly(x1, degree = c1)
    b2 = poly(x2, degree = c2)
    b3 = poly(x1*x2, degree = c3)
    res = cbind(1, b1, b2, b3)
  }
  return(res)
}


series_df_cust = function(X,Y,x_pred,df, type = "ls", std=FALSE, dummy_y =NULL){
  
  x_train = X
  x_pred = x_pred  
  
  y_train = Y
  inv = ginv(t(x_train)%*%x_train)
  if(std==TRUE){
    std = sqrt(x_pred%*%inv %*%t(x_pred ))*sd(Y)
    if (type == "ls"){
      coef = inv %*% t(x_train) %*% y_train
      return (list(y_pred = x_pred %*% coef, sd))
    }else if(type == "forster"){
      sherman_inv = function(x) inv - inv%*% x%*%t(x)%*%inv/ as.numeric(1+ t(x)%*%inv%*%x)
      weight_hn = apply(x_pred, 1, function(x) t(x)%*%sherman_inv(x)%*%x )
      latter = apply(x_pred, 1, function(x) t(x)%*%sherman_inv(x)%*%t(x_train) %*% y_train)
      return (list(as.numeric( 1- weight_hn ) * latter, std))
    }else{print("Type not supported--please input either ls or forster estimator!")}
  }else{
    if (type == "ls"){
      coef = inv %*% t(x_train) %*% y_train
      return (list(y_pred = x_pred %*% coef))
    }else if(type == "forster"){
      sherman_inv = function(x) inv - inv%*% x%*%t(x)%*%inv/ as.numeric(1+ t(x)%*%inv%*%x)
      weight_hn = apply(x_pred, 1, function(x) t(x)%*%sherman_inv(x)%*%x )
      latter = apply(x_pred, 1, function(x) t(x)%*%sherman_inv(x)%*%t(x_train) %*% y_train)
      return (list(as.numeric( 1- weight_hn ) * latter))
    }else{print("Type not supported--please input either ls or forster estimator!")}
  }
}

series_cv_cust = function(X,Y, type = "forster", basis_type = "bs", df_grid = NULL){
  n = dim(X)[1] 
  s_test = floor(n/log(n))
  
  if (!is.null(df_grid))
    df_grid = df_grid
  else if(basis_type == "bs"){
    df_grid = seq(3,15, by = 2)
  }
  else{
    df_grid = seq(1,9, by = 2)
  }
  
  df_grid = expand.grid(df_grid,df_grid,df_grid)
  mse_grid_cv = rep(NA,dim(df_grid)[1])
  
  for (i in seq(dim(df_grid)[1])){
    # cat(paste0(k,"/", KK,' '))
    ind_test = sample(1:n, size = s_test) 
    
    g(train_X, train_Y) %=% list(X[-ind_test,],  Y[-ind_test])
    g(test_X, test_Y) %=% list(X[ind_test,],  Y[ind_test])
    
    basis_train_X = create_basis1(train_X[,1],train_X[,2], k = df_grid[i,], basis_type = basis_type)
    basis_test_X = create_basis1(test_X[,1],test_X[,2], k = df_grid[i,], basis_type = basis_type)
    
    temp = series_df_cust(basis_train_X, train_Y, basis_test_X, type = "forster")[[1]]
    mse_grid_cv[i] = mean((temp - test_Y)^2)
  }  
  ind = which.min(mse_grid_cv)
  
  final_df_grid = cbind(df_grid, mse_grid_cv)
  final_df_grid = final_df_grid[order(mse_grid_cv),]
  return (list(final_df_grid, ind))
}

# ---------

# get a single run of the specified sample size

single_gen = function(n){
  # both x1 and x2 are used, data generation mechanism is more complicated
  
  x1 <- runif(n,0,1);
  x2 <- runif(n,0,1);
  u <- rnorm(n,4,0.5);
  
  # P(Z=1|X)
  prob_z = expit(-1+x1+x2)
  
  z <- rbinom(n,1,prob= prob_z)
  
  # ps_z0xu = exp(alpha_2(U,X))
  ps_z0xu <- exp(-x1 - x2 - u/4)
  
  
  # ps_z1xu = exp(alpha_z(X) + alpha_2(U,X))
  ps_z1xu <- exp(0.5 + (x1 + x2)/2) *  ps_z0xu
  
  a = c(rbinom(n,1,ps_z1xu*z+ps_z0xu*(1-z)))
  
  # E[Y^1|U,Z,X] = E[Y^1|U,A=1,Z,X] = E[Y|U,A=1,Z,X]
  e_y1 =  (x1 + x2 + x1*x2 + z) * exp(u/4)
  eps1 = rnorm(n,0,0.5)
  y1 = e_y1 + eps1
  
  # E[Y^0|U,Z,X] = E[Y^0|U,A=0,X] = E[Y|U,A=0,X]
  e_y0 = (x1 + x2) * exp(u/6)
  eps0 = rnorm(n,0,0.5)
  y0 = e_y0 + eps0
  
  y = a*y1 + (1-a)*y0
  eps = a*eps1 + (1-a)*eps0
  df_pool = data.frame(y,a,z,x1,x2,u,y1,y0,eps,e_y1,e_y0)
  
  return(df_pool)
}


expit <- function(x) {
  exp(x) / (1 + exp(x))
}



### QIV Modifications
single_gen_QIV = function(n){
  ###
  # Simulation Parameters 2/26
  # 
  ###
  
  # Covariates:
  x1 <- runif(n,0,1);
  x2 <- runif(n,0,1);
  u <- rnorm(n,4,0.5);
  
  # P(Z=1|X)
  prob_z = expit(-1+x1+x2)
  
  z <- rbinom(n,1,prob= prob_z)
  
 
  a_func = function(z,u,x1,x2){
    # ps_z0xu = exp(alpha_2(U,X))
    ps_z0xu <- exp(-x1 - x2 - u/4)
    
    # ps_z1xu = exp(alpha_z(X) + alpha_2(U,X))
    ps_z1xu <- exp(0.5 + (x1 + x2)/2) *  ps_z0xu
    pa_x = ps_z1xu*z+ps_z0xu*(1-z)
    return(list("pa_x" = pa_x,"ps_z0xu"=ps_z0xu,"ps_z1xu"=ps_z1xu))
  }
  gen.a = a_func(z,u,x1,x2)
  pa_x=gen.a[["pa_x"]]
  ps_z1xu=gen.a[["ps_z1xu"]]
  ps_z0xu=gen.a[["ps_z0xu"]]
  
  a = c(rbinom(n,1,pa_x))
  
  # Outcome model now relaxes exclusion restriction
  y_func = function(a,z,u,x1,x2){
    beta_z = 2
    # E[Y^1|U,Z,X] = E[Y^1|U,A=1,Z,X] = E[Y|U,A=1,Z,X]
    y1 = function(z,u,x1,x2){
      return((x1 + x2 + x1*x2 ) * exp(u/4) + 2*z)
    }
    e_y1 =  y1(z,u,x1,x2)
    
    # E[Y^0|U,Z,X] = E[Y^0|U,A=0,X] = E[Y|U,A=0,X]
    y0 = function(z,u,x1,x2){
      return((x1 + x2) * exp(u/6) + 2*z)
    }
    e_y0 =  y0(z,u,x1,x2)
  
    phi_1 = y1(z=1,u,x1,x2) - y1(z=0,u,x1,x2)
    ### E(Y|Z=1) = E(Y|Z=1,A=1)*P(A=1|Z=1) + E(Y|Z=1,A=0)*P(A=0|Z=1)
    eyz1 = y1(z=1,u,x1,x2)*ps_z1xu + y0(z=1,u,x1,x2)*(1-ps_z1xu)
    eyz0 = y1(z=0,u,x1,x2)*ps_z0xu + y0(z=0,u,x1,x2)*(1-ps_z0xu)
    delta_y = eyz1 - eyz0 
    y = a*e_y1 + (1-a)*e_y0 
    return(list("y"=y, "e_y1"=e_y1, "e_y0"=e_y0,"phi_1"=phi_1,"delta_y"=delta_y))
  }
  
  
  #df_pool = data.frame(y,a,z,x1,x2,u,y1,y0,eps,e_y1,e_y0)
  
  ### Testing
  
  ### Exclusion Restriction Violation E(Y|A=1,Z=1)-E(Y|A=1,Z=0)
  
  
  gen_y = y_func(a,z,u,x1,x2)   
  e_y0 = gen_y[["e_y0"]]
  e_y1 = gen_y[["e_y1"]]
  eps0 = rnorm(n,0,0.5)
  eps1 = rnorm(n,0,0.5)
  phi_1 = gen_y[["phi_1"]]
  y0 = e_y0 + eps0 
  y1 = e_y1 + eps1
  delta_a = ps_z1xu - ps_z0xu
  delta_y = gen_y[["delta_y"]]
  psi = (delta_y - phi_1)/delta_a
  eps = a*eps1 + (1-a)*eps0
  y = a*y1 + (1-a)*y0
  # Conditional ETT
  ###
  df <- data.frame("y"=y,
                   "a"=a,
                   "z"=z,
                   "x1"=x1,
                   "x2"=x2,
                   "u"=u,
                   "y1"=y1,
                   "y0"=y0,
                   "eps"=eps,
                   "e_Y1"=e_y1,
                   "e_Y0"=e_y0,
                   "phi_1"=phi_1,
                   "psi" = psi,
                   "pi_1" = ps_z1xu,
                   "pi_0" = ps_z0xu, 
                   "pa_x"= pa_x,
                   "pz_x"= prob_z)
  return(df)
}

### Test Assumptions (in case)
test_gen = function(a,z,u,x1,x2,y_func,n){
  
  # assumption 6: exclusion restriction
  # violation of Z on Y 
  # does not depend on A or U 
  
  # (a,z,u,x1,x2)

  y.a0z1i = y_func(0,1,u,x1,x2)[["y"]] 
  y.a0z0i = y_func(0,0,u,x1,x2)[["y"]]
  u.alt = sample(u)
  y.a0z1j = y_func(0,1,u.alt,x1,x2)[["y"]] 
  y.a0z0j = y_func(0,0,u.alt,x1,x2)[["y"]]
  
   
  phi_1.i <- mean(y.a0z1i-y.a0z0i)
  phi_1.j <- mean(y.a0z1j-y.a0z0j)
  tol = 1e-5
  test_6 <- abs(phi_1.i - phi_1.j) <= tol
  if(!(test_6)){
    message = paste("estimated difference in direct effects. Phi_1=",phi_1.i,", ",phi_1.j)
    stop(paste("assumption 6 violated!",message))
  }
  return(y.a0z1i-y.a0z0i)
  
}


### Input: 
###   Y,A,Z,X,
###   Nuissance functions
###   
### From nuissance functions, takes the empirical mean 
### In theory, is centered on true ATT
QIV_IF = function(y,a,z,pi_1,pi_0,pz_x,pa_x,paz_x,phi_1,phi_0,psi,theta){
  p_a = mean(a)
  delta_a = pi_1-pi_0
  fz_x = z*pz_x + (1-z)*(1-pz_x)
  nu_1 = pa_x/p_a*(2*z-1)/fz_x/delta_a
  nu_2 = y - a*psi - z*phi_1 - theta
  nu_3 = pa_x/p_a*a*(2*z-1)/paz_x/delta_a
  nu_4 = y - z*phi_1 - phi_0
  nu_6 = (a*psi/p_a)
  res = (nu_1*nu_2 - nu_3*nu_4 + nu_6)
  
  
  return(res)
}

# Estimate nuisance functions for data_pred
# Returns (conditioned on known covariates X): 
#     "phi_1" = exclusion restriction violation,
#     "phi_0" = baseline,
#     "pi_1"=probability of a=1 given z=1,
#     "pi_0"=probability of a=1 given z=0,
#     "pz_x"=marginal probability of Z=1,
#     "pa_x"=marginal probability of A=1,
#     "paz_x"= joint probability of A=a,Z=z 
QIV_SL <- function(data_train, data_pred,choice = c(1,7)){
  SL.hpara <- list()
  SL.hpara$SLL <- choice
  # Superlearner basic learning algorithms:
  # 1: GLM
  # 2: lasso/ridge
  # 3: earth
  # 4: GAM
  # 5: xgboost
  # 6: polynomial spline
  # 7: random forest
  # 9: gbm
  # 10: 1-layer MLP
  SL.hpara$MLPL <- c(1)
  SL.hpara$MTRY <- c(1)
  SL.hpara$NMN <- 50
  SL.hpara$MLPdecay <- 10^c(-1,-3)
  ###
  X.pos <- which( substr(colnames(data_train),1,1) == "x" )
  Z.pos <- which( substr(colnames(data_train),1,1) == "z" )
  A.pos <- which( substr(colnames(data_train),1,1) == "a" )
  Y.pos <- which( substr(colnames(data_train),1,1) == "y" )
  
  data_train = data.frame(data_train)
  data_pred = data.frame(data_pred)
  
  model_y_za = MySL(Data = data_train,    # training dataset
                    locY = Y.pos,            # response variable Y
                    locX = c(X.pos,Z.pos,A.pos),       # explanatory variables X,Z,A
                    Ydist = gaussian(),      # Y is continuous
                    SL.list = SL.hpara$SLL,  # Machine learning algorithms
                    MTRY = SL.hpara$MTRY,        
                    MLPL = SL.hpara$MLPL,
                    NMN = SL.hpara$NMN,
                    MLPdecay = SL.hpara$MLPdecay)
  ###
  
  Ey_a1z1 = predict(model_y_za, data.frame(data_pred[,X.pos],"z"=1,"a"=1))[[1]]
  Z.model = glm(z ~ .,family = binomial(),data=data_train[,c(Z.pos,X.pos)])
  a_Z.model = glm(a ~ .,family = binomial(),data=data_train[,c(A.pos,Z.pos,X.pos)])
  
  phi_0 = predict(model_y_za, data.frame(data_pred[,X.pos],"z"=0,"a"=1))[[1]]
  phi_1 = Ey_a1z1-phi_0
  Y_z=data_pred[,Y.pos]- data_pred[,Z.pos]*phi_1
  ### Model 1
  #print("Moment Restriction")
  CF <- 2       # cross fitting fold; just use 2 for computational efficiency
  NumCV <- 5    # number of cross validation fold; 5 or 10
  NumCVRep <- 1 # repetition of CV; just use 1 for computational efficiency 
  df_cmr=CMR(Y=Y_z,X = data_pred[,X.pos],Z = data_pred[,Z.pos],A = data_pred[,A.pos],CF,NumCV,NumCVRep)
  psi = df_cmr[["psi"]]
  theta = df_cmr[["theta"]]
  ### Model 2 
  A.p = data_pred[,A.pos]
  Z.p = data_pred[,Z.pos] 
  
  z1.pred = data.frame("z"=1,data_pred[,X.pos])
  pi_1 = predict(a_Z.model,z1.pred,type = c("response"))
  
  z0.pred = data.frame("z"=0,data_pred[,X.pos])
  pi_0 = predict(a_Z.model,z0.pred,type = c("response"))
  
  pz_x = predict(Z.model,data_pred[,X.pos],type = c("response"))
  
  ### marginalize over z 
  # P(A=1|L) = P(Z=1|L)P(A=1|Z=1,L) + P(Z=0|L)P(A=1|Z=0,L)
  pa_x = pz_x*pi_1 + (1-pz_x)*pi_0
  
  
  ### P(A,Z|L) = P(A|Z,L)/P(Z|L)
  p11 = pi_1/pz_x
  p10 = pi_0/(1-pz_x)
  p01 = (1-pi_1)/pz_x
  p00 = (1-pi_0)/(1-pz_x)
  
  paz_x = A.p*Z.p*p11 + A.p*(1-Z.p)*p10 + (1-A.p)*Z.p*p01 + (1-A.p)*(1-Z.p)*p00
  
  
  
  return(list("phi_1" = phi_1,
              "phi_0" = phi_0,
              "pi_1"=pi_1,
              "pi_0"=pi_0,
              "pz_x"=pz_x,
              "pa_x"=pa_x,
              "paz_x"=paz_x,
              "psi"=psi,
              "theta" = theta))
}


CMR = function(Y,X,Z,A,CF,NumCV,NumCVRep){
  ## PMMR regularization parameter
  PL <- -8;      PU <- 0
  ## PMMR bandwidth parameter for the target function
  BW.L <- 0;     BW.U <- 3
  
  Para.Grid <- expand.grid(0, 
                           2,
                           seq(PL,PU,by=2)) 
  ################################################################################
  # Cross-fitting and Cross Validation Sets
  ################################################################################
  N = length(Y)
  
  SS.Index <- list()
  for(ss in 1:CF){
    SS.Index[[ss]] <- ((seq(0,N,length=CF+1))[ss]+1):(seq(0,N,length=CF+1)[ss+1])
  } 
  CV.Index <- list()
  CV.CUT <- seq(0,round(N/CF),length=NumCV+1)
  for(cv in 1:NumCV){
    CV.Index[[cv]] <- ((CV.CUT)[cv]+1):(CV.CUT[cv+1])
  }
  
  if(NumCVRep>1){
    REORDER <- NULL
    for(cv in 1:NumCV){
      REORDER <- rbind(REORDER,((CV.CUT)[cv]+1):(CV.CUT[cv+1]))
    }
    for(tt in 2:NumCVRep){
      SFF <- sapply(1:(N/2/NumCV),function(t){sample(5,5)})
      for(cv in 1:NumCV){
        CV.Index[[cv + (tt-1)*NumCV]] <- 
          sapply(1:(N/2/NumCV),function(t){REORDER[SFF[cv,t],t]})
      }  
    }
  }
  AX <- cbind(A,X)
  ZX <- cbind(Z,X)
  
  Y.MM  <- list()
  A.MM  <- list()
  X.MM  <- list()
  Z.MM  <- list()
  AX.MM <- list()
  ZX.MM <- list()
  
  ### Cross fitting and Validation
  
  for(ss in 1:CF){
    Y.MM [[ss]]  <- Y[SS.Index[[ss]] ]
    A.MM [[ss]]  <- A[SS.Index[[ss]] ]
    X.MM [[ss]]  <- X[SS.Index[[ss]], ]
    Z.MM [[ss]]  <- Z[SS.Index[[ss]] ]
    AX.MM[[ss]]  <- cbind(A.MM[[ss]],X.MM[[ss]])
    ZX.MM[[ss]]  <- cbind(Z.MM[[ss]],X.MM[[ss]])
  }
  Opt.Para.h <- list()
  Opt.Para.h[[1]] <- Opt.Para.h[[2]] <- c(0,1,-2)
  print("CMR")
  h.MM <- list()
  h.predict <- list()
  #browser()
  for(ss in 1:CF){
    h.MM[[ss]] <- 
      FT_PMMR( Y         =Y.MM[[ss]],
               Perturb   =ZX.MM[[ss]],
               Target    =AX.MM[[ss]],
               Diagonal  =rep(1,length(Y.MM[[ss]])),
               Perturb.bw=exp(Opt.Para.h[[ss]][1]),
               Target.bw =exp(Opt.Para.h[[ss]][2]),
               lambda    =exp(Opt.Para.h[[ss]][3]),
               NV        =FALSE)
  }
  
  h.predict[[1]] <- function(AX.New.Input){
    FT_RBF(X     = AX.MM[[1]],
           X.new = AX.New.Input,
           bw.median = exp(Opt.Para.h[[1]][2]))%*%h.MM[[1]]$alpha + h.MM[[1]]$intercept
  }
  
  h.predict[[2]] <- function(AX.New.Input){
    FT_RBF(X     = AX.MM[[2]],
           X.new = AX.New.Input,
           bw.median = exp(Opt.Para.h[[2]][2]))%*%h.MM[[2]]$alpha + h.MM[[2]]$intercept
  }
  
  A1X.MM <- A0X.MM <- AX.MM
  for(ss in 1:CF){
    A1X.MM[[ss]][,1] <- 1
    A0X.MM[[ss]][,1] <- 0
  }
  
  h.A1.hat.CF <- h.A0.hat.CF <- rep(0,N)
  for(ss in 1:CF){
    h.A1.hat.CF[SS.Index[[3-ss]]] <- h.predict[[ss]](A1X.MM[[3-ss]])
    h.A0.hat.CF[SS.Index[[3-ss]]] <- h.predict[[ss]](A0X.MM[[3-ss]])
  }
  
  psi <- h.A1.hat.CF-h.A0.hat.CF
  theta <- h.A0.hat.CF
  return(list("psi"=psi, "theta"=theta))
}

### Input: Y,A,Z,X
### Output: Wald ratio estimator for the ATT 
### Uses nuisance functions obtained by the superlearner
### Constructs a Wald ratio from (this is just psi(L))
QIV_wald = function(a,psi){
  ### E(A * psi(L))/P(A==1)
  return(mean(psi*a)/mean(a==1))
  
}


### Testing
three_fold_eval_QIV = function(data_input, CV = F, k=4, bs = 'bs', choice = c(1,7)){
  y = data_input$y
  a = data_input$a
  z = data_input$z
  x1 = data_input$x1
  x2 = data_input$x2
  s = data_input$s
  
  Data = data.frame(x1,x2,a,z,y,s)
  
  
  # Helper function. Folds have 2 trains and 1 pred. 
  # Evaluates EIF, one_step, and wald on pred.
  # Calls:  
  #   QIV_SL: uses superlearner to estimate nuisance functions
  #   QIV_IF
  #   QIV_wald: constructs wald estimate from psi and a
  
  get_one = function(combination){
    combination = c(1,2,3)
    # choice = c(1,7)
    # CV = T
    # bs = 'bs'
    # k = 4
    
    train_1 = combination[1]
    train_2 = combination[2]
    pred = combination[3]
    
    nuisance_df = QIV_SL(data_train = Data[s==train_1 | s==train_2,],data_pred = Data[s==pred,], choice = choice)

    estimate_wald = a[s==pred]*nuisance_df[["psi"]]/mean(a[s==pred])
    
    
    #y,a,z,pi_1,pi_0,pz_x,pa_x,paz_x,phi_1,phi_0,psi,theta
    
  
    
    estimate_IF = QIV_IF(y = y[s==pred],
                         a = a[s==pred],
                         z = z[s==pred],
                         pi_1 = nuisance_df[["pi_1"]],
                         pi_0 = nuisance_df[["pi_0"]],
                         pz_x = nuisance_df[["pz_x"]],
                         pa_x = nuisance_df[["pa_x"]],
                         paz_x = nuisance_df[["paz_x"]],
                         phi_1 = nuisance_df[["phi_1"]],
                         phi_0 = nuisance_df[["phi_0"]],
                         psi = nuisance_df[["psi"]],
                         theta = nuisance_df[["theta"]])
    
    ### estimates for the ATT (scalar)
    return(data.frame("wald"=estimate_wald,
                      "one_step"=estimate_IF))
    
  }
  
  # s=1,2 used as the train data, s=3, used as the prediction data
  pred_1 = get_one(c(1,2,3))
  # s=1,3 used as the train data, s=2, used as the prediction data
  pred_2 = get_one(c(3,1,2))
  # s=2,3 used as the train data, s=1, used as the prediction data
  pred_3 = get_one(c(2,3,1))
  
  pred_res = rbind(pred_3,pred_2,pred_1)
  return(pred_res)
}

