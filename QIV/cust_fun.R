source("cust_MySL.R")
library(splines)
library(MASS)
library(mgcv)

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
IF = function(y,a,z,p1,p0,e1,e0){
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

three_fold_eval = function(data_input, CV = F, k=4, bs = 'bs', choice = c(1,7)){
  y = data_input$y
  a = data_input$a
  z = data_input$z
  x1 = data_input$x1
  x2 = data_input$x2
  s = data_input$s
  
  Data = data.frame(x1,x2,a,z,y,y*(1-a),s)
  
  get_one = function(combination){
    combination = c(1,2,3)
    # choice = c(1,7)
    # CV = T
    # bs = 'bs'
    # k = 4
    
    train_1 = combination[1]
    train_2 = combination[2]
    pred = combination[3]
    
    nuisance_df_one_step = my_SL(data_train = Data[s==train_1 | s==train_2,],
                                 data_pred = Data[s==pred,], choice = choice)
    nuisance_df_two_step = my_fw(data_train = Data[s==train_1 | s==train_2,],
                                 data_pred = Data[s==pred,], CV = CV, k = k, bs=bs, choice = choice)
    
    estimate_wald = with(nuisance_df_one_step,  (y[s==pred] +  (e1-e0)/(p1-p0))
                         * a[s==pred] ) / mean(a[s==pred])
    estimate_IF = with(nuisance_df_one_step, IF(y[s==pred],a[s==pred],z[s==pred],p1,p0,e1,e0)) / mean(a[s==pred])
    estimate_IF_fw = with(nuisance_df_two_step, IF_fw(y[s==pred],a[s==pred],z[s==pred],
                                                      p1,p0,e1,e0,delta_a_inv_pred, delta_pred)) / mean(a[s==pred])
    estimate_IF_fw2 = with(nuisance_df_two_step, IF_fw2(y[s==pred],a[s==pred],z[s==pred],
                                                        p1,p0,e1,e0,delta_a_inv_pred, delta_pred)) / mean(a[s==pred])
    
    # estimate_IF_fw = 0
    return(data.frame(estimate_wald,
                      estimate_IF,
                      estimate_IF_fw,
                      estimate_IF_fw2))
    
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


