'%=%' = function(l, r, ...) UseMethod('%=%')

#' %=%.lbunch
#' 
#' Given two Grouping objects, copies each item from the list of
#' r into the list in l 
#' 
#' @param l: the Grouping object to be set 
#' @param r: the Grouping object to get from 
#' @param ... 
#'
#' @return none
#' @export
#'
#' @examples
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


#' extendToMatch
#' 
#' Used if LHS(destin) is larger than RHS(source). If the destin argument
#' is a length 1 Grouping object, assume that the only item in the 
#' Grouping object's list is its length. The function returns a 
#' repeated version of source of length equal to length(destin) or destin
#' 
#' @param source: a Grouping object to get from
#' @param destin: a Grouping object to set 
#'
#' @return source: a altered version of the source Grouping object 
#'                 such that the lengths of the underlying lists match
#'
#' @export
#'
#' @examples
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

#' g
#'
#' Constructor for the Grouping/lbunch class used throughout. 
#' Assigns to List a list containing all the function arguments
#' except for the first one.
#'
#' @param ... 
#'
#' @return List
#' @export
#'
#' @examples
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
  
}


#' series_df
#' 
#' Takes in observed covariates, observed outcomes, and
#' kernel basis + degrees of freedom. Returns the predicted Y 
#' on the missing covariates using the regression supplied as 'type'.
#' FW-learner algorithm is implemented with the sherman inverse.
#' type: 
#'    'ls' = least squares, 
#'    'forster' = FW-learner 
#' ("type=forster") basis_type:
#'    'ns' = natural cubic spline of df = df,
#'    'bs' = B-spline of df = df,
#'    'poly' = orthogonal polynomial basis of degree = df
#'    
#'
#' @param X: (n,1) vector of covariates for full data (L1)
#' @param Y: (n,1) vector of responses  
#' @param x_pred: (m,1) vector of covariates for missing data (L2)  
#' @param df: the degrees of freedom for the natural cubic spline basis,
#'            Note: knots = k = df − 1
#' @param type: method of optimization. 
#' @param basis_type: basis function defining the vector space optimized over
#' @param std: if true, returns the standard deviation = 
#'            sqrt(x_pred%*%inv %*%t(x_pred ))*sd(Y) 
#'
#' @return the (m,1) predicted response vector Y^* on x_pred. 
#' @export
#'
#' @examples
series_df = function(X,Y,x_pred,df, type = "ls", basis_type = "poly", std=FALSE){
  or = order(X)
  X = X[or]
  Y = Y[or]
  ll = length(x_pred)
  x_pred = sort(x_pred)
  
  
  if (basis_type == "ns") {
    x_train = ns(X, df )
    x_pred[ll+1] = min(X)
    x_pred[ll+2] = max(X)
    x_pred = matrix(ns(x_pred, df )[1:ll,], nrow = ll)
    
  } else if (basis_type == "poly") {
    x_train = poly(X, df,raw = TRUE)
    x_pred = matrix(poly(x_pred, df,raw = TRUE),nrow = ll)
  } else if (basis_type == "bs") {
    x_train = bs(X, df = 4, knots = df)
    x_pred[ll+1] = min(X)
    x_pred[ll+2] = max(X)
    x_pred = matrix(bs(x_pred, df = 4, knots = df)[1:ll,], nrow = ll)
  } else {
    print("Basis type not supported: please use either bs, ns or poly")
  }
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

#' series_cv_new
#' 
#' Method for finding optimal degrees of freedom argument(df) 
#' for the FW learner for a given set of covariates, responses,
#' missing data, and basis_type. Uses cross validation on KK 
#' iterations on n/log(n) sample points from n
#' 
#'
#' @param X: (vector) of covariates
#' @param Y: (vector) vector of responses
#' @param x_pred: (vector) of covariates for missing data 
#' @param type: (string) type supplied to series_df 
#' @param basis_type: (string) basis_type supplies to series_df
#' @param KK: (Int) Number of iterations to run cross validation
#' @param df: (Int) "df" supplied to series_df 
#' @param std: (Boolean) "std" supplied to series_df 
#' @param df_grid: (vector) of df values to check. 
#'
#' @return List(final_est, estimator_list,D_list): 
#'         final_est: the mean estimator_list
#'         estimator_list: the predicted values from series_df for 
#'                         each split with the minimum MSE across all 
#'                         df in df_grid  
#'         D_list: list of the df in df_grid corresponding to the predicted
#'                 values in estimator_list for split k < KK
#' @export
#'
#' @examples
series_cv_new = function(X,Y,x_pred, type = "forster", basis_type = "poly", KK=1, df = NULL, std=FALSE, df_grid = NULL){
  if (is.null(std)){std=FALSE}
  if (!(is.null(df))){
    return (series_df(X,Y,x_pred, df = df, type = type, basis_type = basis_type, std=std))
  }
  n = length(X) 
  s_test = floor(n/log(n))
  if (is.null(df_grid)){
    df_grid = seq(4,18, by = 2)
    if (basis_type == "bs"){df_grid = seq(50,200, by = 15)}
  }
  l = length(df_grid)
  mse_grid_cv = matrix(NA,l, KK)
  estimator_list = D_list = rep(NA,KK)
  for (k in seq(KK)){
    cat(paste0(k,"/", KK,' '))
    ind_test = sample(1:n, size = s_test) 
    
    g(train_X, train_Y) %=% list(X[-ind_test],  Y[-ind_test])
    g(test_X, test_Y) %=% list(X[ind_test],  Y[ind_test])
    for (i in 1:l){
      df_temp = df_grid[i]
      temp =  series_df(X = train_X, Y=train_Y, x_pred = test_X, df = df_temp ,type = type, basis_type = basis_type,std=std)[[1]]
      mse_grid_cv[i,k] = mean((temp - test_Y)^2)
    }
    ind = which(mse_grid_cv[,k]==min(mse_grid_cv[,k]))[1]
    D_list[k] = df_grid[ind]
    estimator_list[k] =  series_df(X = X, Y=Y, x_pred = x_pred, df = D_list[k] ,type = type, basis_type = basis_type,std=std)[[1]]
  }  
  
  return (list(final_est = mean(estimator_list), estimator_list = estimator_list, D_list = D_list))
}

# nonparametric RKHS fitted DR estimator for CATE through Cross-fitting
#' CF_cate
#'  
#' Used to generate pseudo-outcomes for the proximal causal inference setting 
#' 
#' @param data:  dataframe containing columns X,Z,W,A,Y
#' @param kernel_sigma_h: sigma parameter for the radial basis kernel for the h functions
#' @param kernel_sigma_q: sigma parameter for the radial basis kernel for the q functions
#' @param lm_Hh1: 
#' @param lm_Qh1:  
#' @param lm_Hh0:  
#' @param lm_Qh0:  
#' @param lm_Hq1:  
#' @param lm_Qq1:  
#' @param lm_Hq0:  
#' @param lm_Qq0:  
#' @param CF_K_fold: number of folds to split data into for cross validation 
#' @param eval_data: if supplied, then use this as the evaluation(full data) set for cross validation
#' @param center 
#'
#' @return list(DR_IF, DR_est, sqrt(sum(DR_IF^2)))
#' @export
#'
#' @examples
CF_cate <- function(data, kernel_sigma_h, kernel_sigma_q, lm_Hh1, lm_Qh1, lm_Hh0, lm_Qh0, lm_Hq1, lm_Qq1, lm_Hq0, lm_Qq0, CF_K_fold, eval_data=NULL, center=FALSE){
  n_sample <- length(data$Y)
  working_data <- with(data, list(X = X, Z = Z, W = W, A = A, Y = Y))
  random_sample <- createFolds(1:n_sample, k = CF_K_fold)
  OR_est <- 0
  IPW_est <- 0
  DR_est <- 0
  OR_IF <- c()
  IPW_IF <- c()
  DR_IF <- c()
  if (is.null(eval_data)){
    for (cf_rep in 1:CF_K_fold) {
      tr_data <- sapply(working_data, function(x) matrix(x[-random_sample[[cf_rep]], ], ncol = ncol(x)))
      eval_data <- sapply(working_data, function(x) matrix(x[random_sample[[cf_rep]], ], ncol = ncol(x)))
      n_tr <- n_sample - length(random_sample[[cf_rep]])
      n_eval <- length(random_sample[[cf_rep]])
      
      ## Gaussian RBF kernel ############################################################
      h_arg <- cbind(tr_data$W, tr_data$X)
      q_arg <- cbind(tr_data$Z, tr_data$X)
      K_H <- rbfkernel(h_arg, sigma = kernel_sigma_h)
      K_Q <- rbfkernel(q_arg, sigma = kernel_sigma_q)
      ###################################################################################
      
      ## Optimization ###################################################################
      ident_n <- diag(rep(1, n_tr))
      D1 <- diag(as.numeric(tr_data$A), n_tr, n_tr)
      D0 <- diag(as.numeric(1 - tr_data$A), n_tr, n_tr)
      ## Optimization h ##
      Gamm <- 0.25 * K_Q %*% ginv(K_Q / n_tr + lm_Qh1 / n_tr^0.8 * ident_n)
      alpha1h <- ginv(K_H %*% D1 %*% Gamm %*% D1 %*% K_H + n_tr^2 * lm_Hh1 / n_tr^0.8 * K_H) %*%
        K_H %*% D1 %*% Gamm %*% (tr_data$Y * tr_data$A)
      Gamm <- 0.25 * K_Q %*% ginv(K_Q / n_tr + lm_Qh0 / n_tr^0.8 * ident_n)
      
      alpha0h <- ginv(K_H %*% D0 %*% Gamm %*% D0 %*% K_H + n_tr^2 * lm_Hh0 / n_tr^0.8 * K_H) %*%
        K_H %*% D0 %*% Gamm %*% (tr_data$Y * (1 - tr_data$A))
      ## Optimization q ##
      Gamm <- 0.25 * K_H %*% ginv( K_H / n_tr + lm_Hq1 / n_tr^0.8 * ident_n )
      alpha1q <- ginv(K_Q %*% D1 %*% Gamm %*% D1 %*% K_Q + n_tr^2 * lm_Qq1 / n_tr^0.8 * K_Q) %*%
        K_Q %*% D1 %*% Gamm %*% matrix(1, n_tr, 1)
      #Gamm <- 0.25 * K_H %*% ginv( K_H / n_tr + lm_Hq0 / n_tr^0.8 * ident_n )
      alpha0q <- ginv(K_Q %*% D0 %*% Gamm %*% D0 %*% K_Q + n_tr^2 * lm_Qq0 / n_tr^0.8 * K_Q) %*%
        K_Q %*% D0 %*% Gamm %*% matrix(1, n_tr, 1)
      ###################################################################################
      
      ## Estimate function h_hat ########################################################
      h_hat_1 <- function(w1, x1){
        T_mat <- c(w1,x1)
        S <- rbfkernel(h_arg, sigma = kernel_sigma_h, matrix(T_mat,1,length(T_mat)))
        return(t(S) %*% alpha1h)
      }
      h_hat_0 <- function(w1, x1){
        T_mat <- c(w1,x1)
        S <- rbfkernel(h_arg, sigma = kernel_sigma_h, matrix(T_mat,1,length(T_mat)))
        return(t(S) %*% alpha0h)
      }
      h_hat <- function(w1, a1, x1){
        return(a1*h_hat_1(w1,x1)+(1-a1)*h_hat_0(w1,x1))
      }
      ###################################################################################
      
      ## Estimate function q_hat ########################################################
      q_hat_1 <- function(z1, x1){
        T_mat <- c(z1,x1)
        S <- rbfkernel(q_arg, sigma = kernel_sigma_q,matrix(T_mat, 1, length(T_mat)))
        return(t(S) %*% alpha1q)
      }
      q_hat_0 <- function(z1, x1){
        T_mat <- c(z1,x1)
        S <- rbfkernel(q_arg, sigma = kernel_sigma_q,matrix(T_mat, 1, length(T_mat)))
        return(t(S) %*% alpha0q)
      }
      q_hat <- function(z1, a1, x1){
        return(a1 * q_hat_1(z1, x1) + (1 - a1) * q_hat_0(z1, x1))
      }
      ###################################################################################
      ## Evaluation DR ##################################################################
      est <- 0
      var_est <- 0
      for (j in seq(n_eval)){
        if (j %% 100 ==0){cat(j, 'out of', n_eval,';')}
        IF <-   (-1)^(1 - eval_data$A[j, ]) * q_hat(eval_data$Z[j, ], eval_data$A[j, ], eval_data$X[j, ]) *
          (eval_data$Y[j, ] - h_hat(eval_data$W[j, ], eval_data$A[j, ], eval_data$X[j, ])) +
          (h_hat(eval_data$W[j, ], 1, eval_data$X[j, ]) - h_hat(eval_data$W[j, ], 0, eval_data$X[j,]))
        DR_IF <- c(DR_IF, IF)
      }
    }  
    ###################################################################################
    DR_IF <- DR_IF / CF_K_fold
    DR_est <- mean(DR_IF)
    DR_IF <- DR_IF
  }else{
    tr_data <- working_data
    n_tr <- n_sample
    n_eval <- length(eval_data$Y)
    
    ## Gaussian RBF kernel ############################################################
    h_arg <- cbind(tr_data$W, tr_data$X)
    q_arg <- cbind(tr_data$Z, tr_data$X)
    K_H <- rbfkernel(h_arg, sigma = kernel_sigma_h)
    K_Q <- rbfkernel(q_arg, sigma = kernel_sigma_q)
    ###################################################################################
    
    ## Optimization ###################################################################
    ident_n <- diag(rep(1, n_tr))
    D1 <- diag(as.numeric(tr_data$A), n_tr, n_tr)
    D0 <- diag(as.numeric(1 - tr_data$A), n_tr, n_tr)
    ## Optimization h ##
    Gamm <- 0.25 * K_Q %*% ginv(K_Q / n_tr + lm_Qh1 / n_tr^0.8 * ident_n)
    alpha1h <- ginv(K_H %*% D1 %*% Gamm %*% D1 %*% K_H + n_tr^2 * lm_Hh1 / n_tr^0.8 * K_H) %*%
      K_H %*% D1 %*% Gamm %*% (tr_data$Y * tr_data$A)
    Gamm <- 0.25 * K_Q %*% ginv(K_Q / n_tr + lm_Qh0 / n_tr^0.8 * ident_n)
    
    alpha0h <- ginv(K_H %*% D0 %*% Gamm %*% D0 %*% K_H + n_tr^2 * lm_Hh0 / n_tr^0.8 * K_H) %*%
      K_H %*% D0 %*% Gamm %*% (tr_data$Y * (1 - tr_data$A))
    ## Optimization q ##
    Gamm <- 0.25 * K_H %*% ginv( K_H / n_tr + lm_Hq1 / n_tr^0.8 * ident_n )
    alpha1q <- ginv(K_Q %*% D1 %*% Gamm %*% D1 %*% K_Q + n_tr^2 * lm_Qq1 / n_tr^0.8 * K_Q) %*%
      K_Q %*% D1 %*% Gamm %*% matrix(1, n_tr, 1)
    #Gamm <- 0.25 * K_H %*% ginv( K_H / n_tr + lm_Hq0 / n_tr^0.8 * ident_n )
    alpha0q <- ginv(K_Q %*% D0 %*% Gamm %*% D0 %*% K_Q + n_tr^2 * lm_Qq0 / n_tr^0.8 * K_Q) %*%
      K_Q %*% D0 %*% Gamm %*% matrix(1, n_tr, 1)
    ###################################################################################
    
    ## Estimate function h_hat ########################################################
    h_hat_1 <- function(w1, x1){
      T_mat <- c(w1,x1)
      S <- rbfkernel(h_arg, sigma = kernel_sigma_h, matrix(T_mat,1,length(T_mat)))
      return(t(S) %*% alpha1h)
    }
    h_hat_0 <- function(w1, x1){
      T_mat <- c(w1,x1)
      S <- rbfkernel(h_arg, sigma = kernel_sigma_h, matrix(T_mat,1,length(T_mat)))
      return(t(S) %*% alpha0h)
    }
    h_hat <- function(w1, a1, x1){
      return(a1*h_hat_1(w1,x1)+(1-a1)*h_hat_0(w1,x1))
    }
    ###################################################################################
    
    ## Estimate function q_hat ########################################################
    q_hat_1 <- function(z1, x1){
      T_mat <- c(z1,x1)
      S <- rbfkernel(q_arg, sigma = kernel_sigma_q,matrix(T_mat, 1, length(T_mat)))
      return(t(S) %*% alpha1q)
    }
    q_hat_0 <- function(z1, x1){
      T_mat <- c(z1,x1)
      S <- rbfkernel(q_arg, sigma = kernel_sigma_q,matrix(T_mat, 1, length(T_mat)))
      return(t(S) %*% alpha0q)
    }
    q_hat <- function(z1, a1, x1){
      return(a1 * q_hat_1(z1, x1) + (1 - a1) * q_hat_0(z1, x1))
    }
    ###################################################################################
    ## Evaluation DR ##################################################################
    est <- 0
    var_est <- 0
    for (j in seq(n_eval)){
      if (j %% 100 ==0){cat(j, 'out of', n_eval,';')}
      IF <-   (-1)^(1 - eval_data$A[j, ]) * q_hat(eval_data$Z[j, ], eval_data$A[j, ], eval_data$X[j, ]) *
        (eval_data$Y[j, ] - h_hat(eval_data$W[j, ], eval_data$A[j, ], eval_data$X[j, ])) +
        (h_hat(eval_data$W[j, ], 1, eval_data$X[j, ]) - h_hat(eval_data$W[j, ], 0, eval_data$X[j,]))
      DR_IF <- c(DR_IF, IF)
    }
    DR_IF <- DR_IF / CF_K_fold
    DR_est <- mean(DR_IF)
    DR_IF <- DR_IF
    ###################################################################################
  }
  if (center==TRUE){DR_IF = DR_IF - mean(DR_IF)}
  result <- list(DR_IF, DR_est, sqrt(sum(DR_IF^2)))
  return(result)
}
