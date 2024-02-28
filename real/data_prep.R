rm(list = ls())
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
set.seed(8123)
#library(ggplot2)
library(MASS)
library(corpcor)
library(rdetools)
library(expm)
###############################
### DATA ANALYSIS FOR PAPER
############ data process ####################
#library(trust)
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
# source("ate_identity.R")
# source("ate_logit.R")
# source("ate_MLE_identity.R")
# source("ate_MLE_logit.R")
# source("ate_nlm_identity.R")
# source("ate_nlm_logit.R")

expit<-function(x){exp(x)/(1+exp(x))}

### 1: data-preparation
data_prep = function(d = 7){
  data1<-read.table("rhc.csv",sep=",",head=TRUE)
  data1$sex<-ifelse(data1$sex=="Female",1,0) # Female = 1 and Male = 0
  # race        black/white/other, numeric  CATEGORICAL: MAKE BINARY INDICATORS (2)
  levels(data1$race)
  data1$raceblack<-ifelse(data1$race=="black",1,0) # Black = 1 
  data1$raceother<-ifelse(data1$race=="other",1,0) # other = 1    white is reference 
  # edu         years of eduction 
  # income      income categories           CATEGORICAL: MAKE BINARY INDICATORS (3)
  levels(data1$income)
  data1$income1<-ifelse(data1$income=="$11-$25k",1,0)
  data1$income2<-ifelse(data1$income=="$25-$50k",1,0)
  data1$income3<-ifelse(data1$income=="> $50k",1,0)               # under $11k is reference 
  # ninsclas    type of insurance           CATEGORICAL: MAKE BINARY INDICATORS (5)
  levels(data1$ninsclas)
  data1$ins_care<-ifelse(data1$ninsclas=="Medicare",1,0)
  data1$ins_pcare<-ifelse(data1$ninsclas=="Private & Medicare",1,0)
  data1$ins_caid<-ifelse(data1$ninsclas=="Medicaid",1,0)
  data1$ins_no<-ifelse(data1$ninsclas=="No insurance",1,0)
  data1$ins_carecaid<-ifelse(data1$ninsclas=="Medicare & Medicaid",1,0)      # private is reference
  # cat1        primary disease category    CATEGORICAL: MAKE BINARY INDICATORS (8)
  data1$cat1_copd<-ifelse(data1$cat1=="COPD",1,0)
  data1$cat1_mosfsep<-ifelse(data1$cat1=="MOSF w/Sepsis",1,0)
  data1$cat1_mosfmal<-ifelse(data1$cat1=="MOSF w/Malignancy",1,0)
  data1$cat1_chf<-ifelse(data1$cat1=="CHF",1,0)
  data1$cat1_coma<-ifelse(data1$cat1=="Coma",1,0)
  data1$cat1_cirr<-ifelse(data1$cat1=="Cirrhosis",1,0)
  data1$cat1_lung<-ifelse(data1$cat1=="Lung Cancer",1,0)
  data1$cat1_colon<-ifelse(data1$cat1=="Colon Cancer",1,0)    # ARF is reference
  # cat2        secondary disease category  CATEGORICAL: MAKE BINARY INDICATORS (6)      a lot of missingness! extra category? (See Hirano/Imbens)
  data1$cat2NA<-factor(data1$cat2, exclude=NULL) 
  data1$cat2_mosfsep<-ifelse(data1$cat2NA=="MOSF w/Sepsis",1,0)
  data1$cat2_coma<-ifelse(data1$cat2NA=="Coma",1,0)
  data1$cat2_mosfmal<-ifelse(data1$cat2NA=="MOSF w/Malignancy",1,0)
  data1$cat2_lung<-ifelse(data1$cat2NA=="Lung Cancer",1,0)
  data1$cat2_cirr<-ifelse(data1$cat2NA=="Cirrhosis",1,0)
  data1$cat2_colon<-ifelse(data1$cat2NA=="Colon Cancer",1,0)    # NA is reference
  # resp        yes/no                      CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$resp<-ifelse(data1$resp=="Yes",1,0)
  # card        yes/no                      CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$card<-ifelse(data1$card=="Yes",1,0)
  # neuro       yes/no                      CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$neuro<-ifelse(data1$neuro=="Yes",1,0)
  # gastr       yes/no                      CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$gastr<-ifelse(data1$gastr=="Yes",1,0)
  # renal       yes/no                      CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$renal<-ifelse(data1$renal=="Yes",1,0)
  # meta        yes/no                      CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$meta<-ifelse(data1$meta=="Yes",1,0)
  # hema        yes/no                      CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$hema<-ifelse(data1$hema=="Yes",1,0)
  # seps        yes/no                      CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$seps<-ifelse(data1$seps=="Yes",1,0)
  # trauma      yes/no                      CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$trauma<-ifelse(data1$trauma=="Yes",1,0)
  # ortho       yes/no                      CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$ortho<-ifelse(data1$ortho=="Yes",1,0)
  # adld3p      ADL, numeric                                                         a lot of missingness!
  data1$MISSadld3p<-ifelse(is.na(data1$adld3p)=="TRUE",1,0)
  # das2d3pc    Index, numeric         
  # dnr1        yes/no                      CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$dnr1<-ifelse(data1$dnr1=="Yes",1,0)
  # ca          yes/no/meta                 CATEGORICAL: MAKE BINARY INDICATOR (2)
  data1$ca_yes<-ifelse(data1$ca=="Yes",1,0)
  data1$ca_meta<-ifelse(data1$ca=="Metastatic",1,0)
  # surv2md1    estimated surv prob 
  # aps1        score, numeric
  # scoma1      score, numeric
  # wtkilo1     weight in kilo, numeric                                              a lot of zero measurements: make extra indicator saying this is 0
  data1$wt0<-ifelse(data1$wtkilo1==0,1,0)                                                
  data1$MISSurin1<-ifelse(is.na(data1$urin1)=="TRUE",1,0)
  # cardiohx    yes/no, already in 0/1
  # chfhx       yes/no, already in 0/1
  # dementhx    yes/no, already in 0/1
  # psychhx     yes/no, already in 0/1
  # chrpulhx    yes/no, already in 0/1
  # renalhx     yes/no, already in 0/1
  # liverhx     yes/no, already in 0/1
  # gibledhx    yes/no, already in 0/1
  # malighx     yes/no, already in 0/1
  # immunhx     yes/no, already in 0/1
  # transhx     yes/no, already in 0/1
  # amihx       yes/no, already in 0/1
  
  ## Treatment variable
  # swang1      RHC or NO RHC               CATEGORICAL: MAKE BINARY INDICATOR (1)
  data1$swang1<-ifelse(data1$swang1=="RHC",1,0)   # 1 indicates RHC was given, 0 means NO RHC
  table(data1$swang1)
  
  ## Outcome variable
  # dth30       death after 30 days
  data1$dth30<-ifelse(data1$dth30=="No",1,0)    # 1 indicates survival, we want to model the probability of surviving
  table(data1$dth30)
  
  
  Y<-data1$t3d30
  A<-data1$swang1
  attach(data1)
  
  

  
  alpha <- 0.05
  n <- length(Y)
  ############################ data ############################
  W <- cbind(ph1,hema1)
  Z <- cbind(pafi1,paco21) 
  X1 <- cbind(age,sex,raceblack,raceother,edu,income1,income2,income3,ins_care,ins_pcare,ins_caid,ins_no,ins_carecaid,cat1_copd,cat1_mosfsep,
             cat1_mosfmal,cat1_chf,cat1_coma,cat1_cirr,cat1_lung,cat1_colon,cat2_mosfsep,cat2_coma,cat2_mosfmal,cat2_lung,cat2_cirr,
             cat2_colon,resp,card,neuro,gastr,renal,meta,#hema,
             seps,trauma,ortho,das2d3pc,dnr1,ca_yes,ca_meta,surv2md1,aps1,scoma1,wtkilo1,
             temp1,meanbp1,resp1,hrt1,wblc1,sod1,pot1,crea1,bili1,alb1,cardiohx,chfhx,dementhx,psychhx,chrpulhx,renalhx,
             liverhx,gibledhx,malighx,immunhx,transhx,amihx,wt0)
  if (d == 67){
   X = X1
  }else if (d == 71){
    X = cbind(Z,W,X1)
  } else if (d == 7){
    X <- cbind(age, sex, cat1_coma, cat2_coma, dnr1, surv2md1, aps1)  
  }
  data <- list(A=matrix(A), Y=matrix(Y), Z=Z, W=W, X=X)
  attach(data)
}

