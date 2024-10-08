################################################################################
# Required Packages
################################################################################

library(SuperLearner)
library(caret)
library(nnet)
library(glmnet)
library(earth)
library(gam)
library(gbm)
library(xgboost)     
library(kernlab)
library(polspline)
library(ranger)

################################################################################
# MySL : Superlearner-based nuisance function estimation
# 
# This function is used to nonparametrically estimate
# the outcome regression and the propensity score
#
# Input: 
# Data : a dataframe containing independent variables and a dependent variable
# locY : a scalar that indicates the column index of a dependent variable
# locX : a vector that indicates the column indices of independent variables
# Ydist: type of a dependent variable (gaussian()=continuous, binomial()=binary)
# SL.list : Superlearner basic learning algorithms 
#   1: generalized linear model
#   2: lasso/ridge/elastic net regressions
#   3: multivariate adaptive regression splines (earth)
#   4: generalized additive model (GAM)
#   5: xgboost
#   6: polynomial spline
#   7: random forest
#   8 : 1-layer neural net
#   9: gradient boosting method (GBM)
#   10: single-layer perceptron (SLP)
#   11: multi-layer perceptron
# MTRY: hyperparameter for random forest
# MLPL: hyperparameter for 1-layer MLP (number of nodes)
# MLPdecay: hyperparameter for 1-layer MLP (decay parameter)
# NMN: hyperparameter for gbm
# obsWeights: weight vector for each observation
# PS.thr: threshold for truncating the estimated propensity score
# CVlist: number of cross-validation folds
################################################################################

MySL <- function( Data, 
                  locY, 
                  locX, 
                  Ydist=gaussian(), 
                  SL.list=c(1:11), 
                  MTRY=c(2,4,6,8), 
                  MLPL=c(2,4,6,8), 
                  MLPdecay=c(10^(-4),10^(-5)), 
                  NMN=c(20), 
                  obsWeights=NULL,
                  PS.thr = 10^(-3), 
                  CVlist=NULL ){
    
    ## Poisson c(2,4,5)
    if(Ydist$family=="poisson"){
        SL.list <- intersect(c(1,2,4,5),SL.list)
    }
    
    Learners <- list() 
    
    SL.caret.SLP <- function (Y, X, newX, family, obsWeights, method = "mlpML", 
                              L1,L2,L3,decay,
                              trControl = caret::trainControl(method = "none"), 
                              metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...) 
    {
        if (family$family == "gaussian") {
            fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      preProc =  c('center', 'scale', 'pca'),
                                      hiddenActFunc = "Act_Logistic",linOut=TRUE,
                                      learnFuncParams=decay,
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "raw")
        }
        if (family$family == "binomial") {
            Y.f <- as.factor(Y)
            levels(Y.f) <- c("A0", "A1")
            fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      hiddenActFunc = "Act_Identity",
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
        }
        fit <- list(object = fit.train)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.caret")
        return(out)
    }
    
    
    SL.caret.MLP <- function (Y, X, newX, family, obsWeights, method = "mlpML", 
                              L1,L2,L3,decay,
                              trControl = caret::trainControl(method = "none"), 
                              metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...) 
    {
        if (family$family == "gaussian") {
            fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      learnFuncParams=decay,
                                      hiddenActFunc = "Act_Logistic",linOut=TRUE,
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "raw")
        }
        if (family$family == "binomial") {
            Y.f <- as.factor(Y)
            levels(Y.f) <- c("A0", "A1")
            fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      hiddenActFunc = "Act_Identity",
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
        }
        fit <- list(object = fit.train)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.caret")
        return(out)
    }
    
    
    
    SL.caret.gbm <- function (Y, X, newX, family, obsWeights, method = "gbm",
                              ntree,intdepth,sh,nmn,
                              trControl = caret::trainControl(method = "none"),
                              metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...)
    {
        if (family$family == "gaussian") {
            
            fit.train <- caret::train(x = X, y = Y, weights = obsWeights,
                                      metric = metric, method = method,
                                      tuneGrid = expand.grid(n.trees=ntree,interaction.depth=intdepth,shrinkage=sh,n.minobsinnode=nmn),
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "raw")
        }
        if (family$family == "binomial") {
            Y.f <- as.factor(Y)
            levels(Y.f) <- c("A0", "A1")
            fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                      metric = metric, method = method,
                                      tuneGrid = expand.grid(n.trees=ntree,interaction.depth=intdepth,shrinkage=sh,n.minobsinnode=nmn),
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
        }
        fit <- list(object = fit.train)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.caret")
        return(out)
    }
     
    
    SL.new.earth <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                              nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                              nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...)
    {
        if (family$family == "gaussian") {
            fit.earth <- earth::earth(x = X, y = Y, degree = degree, weights=obsWeights, 
                                      nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                                      ncross = ncross, minspan = minspan, endspan = endspan)
        }
        if (family$family == "binomial") {
            fit.earth <- earth::earth(x = X, y = Y, degree = degree, weights=obsWeights,
                                      nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                                      ncross = ncross, minspan = minspan, endspan = endspan, 
                                      glm = list(family = binomial))
        }
        pred <- predict(fit.earth, newdata = newX, type = "response")
        fit <- list(object = fit.earth)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.earth")
        return(out)
    }
    
    SL.new.xgboost <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
                                max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                                nthread = 1, verbose = 0, save_period = NULL, ...) 
    {
        if (packageVersion("xgboost") < 0.6) 
            stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
        if (!is.matrix(X)) {
            X = model.matrix(~. - 1, X)
        }
        xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
        if (family$family == "gaussian") {
            model = xgboost::xgboost(data = xgmat, objective = "reg:linear", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, nthread = nthread, 
                                     params = params, save_period = save_period)
        }
        if (family$family == "binomial") {
            model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, nthread = nthread, 
                                     params = params, save_period = save_period)
        }
        if (family$family == "multinomial") {
            model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                                     nthread = nthread, params = params, save_period = save_period)
        }
        if (family$family == "poisson") {
            model = xgboost::xgboost(data = xgmat, objective = "count:poisson", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, 
                                     nthread = nthread, params = params, save_period = save_period)
        }
        if (!is.matrix(newX)) {
            newX = model.matrix(~. - 1, newX)
        }
        pred = predict(model, newdata = newX)
        fit = list(object = model)
        class(fit) = c("SL.xgboost")
        out = list(pred = pred, fit = fit)
        return(out)
    }
    
    Learners[[1]] <- create.Learner("SL.glm")
    TOTAL.M <- 1
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.glmnet",tune=list(alpha=c(1,0.5,0),useMin=c(TRUE,FALSE))) 
    # Lasso.min , EN.min , Ridge.min , Lasso.1se , EN.1se , Ridge.1se
    TOTAL.M <- TOTAL.M+1           # 2
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.new.earth",tune=list(degree=c(1,2,3,4,5)))
    # Earth.deg=1 , Earth.deg=2 , Earth.deg=3 , Earth.deg=4 , Earth.deg=5
    TOTAL.M <- TOTAL.M+1           # 3
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.gam",tune=list(deg.gam=c(1,2,3,4,5)))
    # Gam.deg=1 , Gam.deg=2 , Gam.deg=3 , Gam.deg=4 , Gam.deg=5
    TOTAL.M <- TOTAL.M+1           # 4
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.new.xgboost",tune=list(n.trees=c(100,300,500),max_depth=c(1,2,3,4)))
    # xgboost
    TOTAL.M <- TOTAL.M+1           # 5
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.polymars",tune=list(knots=c(2,3,4)))
    # polspline (similar to earth)
    TOTAL.M <- TOTAL.M+1           # 6
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.ranger",tune=list(num.trees=c(500,1000,1500),mtry=MTRY))
    # RF
    TOTAL.M <- TOTAL.M+1           # 7
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.nnet",tune=list(linout=c(TRUE,FALSE), decay=c(0,0.1),size=MTRY))
    # nnet
    TOTAL.M <- TOTAL.M+1           # 8
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.caret.gbm",tune=list(ntree=c(100,300,500),intdepth=c(1,2,3),sh=c(0.1,0.01),nmn=NMN))
    # gbm
    TOTAL.M <- TOTAL.M+1           # 9
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.caret.SLP", tune = list(L1=MLPL,L2=c(0),L3=c(0),decay=MLPdecay))
    # 1layer
    TOTAL.M <- TOTAL.M+1           # 10
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.caret.MLP", tune = list(L1=MLPL,L2=MLPL,L3=MLPL,decay=MLPdecay))
    # 3layer
    TOTAL.M <- TOTAL.M+1           # 11

    BaseLearner <- Learners[[ SL.list[1] ]]$names
    
    if( length(SL.list)>1 ){
        for( METHOD in 2:length(SL.list) ){
            BaseLearner <- c(BaseLearner, Learners[[ SL.list[METHOD] ]]$names)
        }
    }
    
    if( length(locX)==1 ){
        dX <- data.frame(matrix( Data[,locX], dim(Data)[1], 1))
        colnames(dX) <- colnames(Data)[locX]
    } else {
        dX <- Data[,locX]
    }
    
    if(is.null(CVlist)){
        CVL <- list(V = 5)
    } else {
        CVL <- list(V = 5, shuffle = FALSE, validRows = CVlist)
    }
    
    
    if(Ydist$family!="binomial"){
        capture.output( Fitted.SL <- SuperLearner(Y=Data[,locY],X=dX,family=Ydist,
                                                  SL.library=BaseLearner, cvControl = CVL, obsWeights=obsWeights) , file=NULL )
        
    } else {
        capture.output( Fitted.SL <- SuperLearner(Y=Data[,locY],X=dX,family=Ydist,
                                                   SL.library=BaseLearner, cvControl = CVL, obsWeights=obsWeights) , file=NULL )
    }
    
    
    return(Fitted.SL)
}