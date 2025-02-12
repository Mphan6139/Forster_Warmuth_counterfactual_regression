################################################################################
# Setup:
# Y = h(A,X)+e
# E[e|A,X] =/=0
# E[e|Z,X] =0
# This is equal to solve E[ (Y-h(A,X))*g(Z,X) ] = 0 for any g
# Target function = f(A,X) with target variable A,X
# Perturbing function = g(Z,X) with perturbing variable (Z,X)
################################################################################

################################################################################
# Basic Simulation Parameters
################################################################################

N <- 2000
CF <- 2       # cross fitting fold; just use 2 for computational efficiency
NumCV <- 5    # number of cross validation fold; 5 or 10
NumCVRep <- 1 # repetition of CV; just use 1 for computational efficiency 

################################################################################
# Package and Source Files
################################################################################

source("0.Functions_PMMR.R")

################################################################################
# Parameters
################################################################################

## PMMR regularization parameter
PL <- -8;      PU <- 0
## PMMR bandwidth parameter for the target function
BW.L <- 0;     BW.U <- 3

Para.Grid <- expand.grid(0, 
                         2,
                         seq(PL,PU,by=2)) 

################################################################################
# DGP
################################################################################

X <- rnorm(N)
U <- rnorm(N) # unmeasured confounder
Z <- rbinom(N,1,expit(0.5*X))
A <- rbinom(N,1,expit(-0.25+0.5*Z+0.25*X+0.25*U))
error <- 0.25*rnorm(N)+0.25*U
Y <- 1+0.25*X+(0.5+X)*A+error

cor(error,A)  ## (incrase N to a large number to see this does not go to zero)
cor(error,Z)  ## (incrase N to a large number to see this goes to zero)

ATT.true <- (0.5+X)

Data <- cbind(Y,A,Z,X)
colnames(Data) <- c("Y","A","Z","X")

# ## make test data for evaluation
# 
# X.Test <- rnorm(N)
# U.Test <- rnorm(N) # unmeasured confounder
# Z.Test <- rbinom(N,1,expit(0.5*X.Test))
# A.Test <- rbinom(N,1,expit(Z.Test+0.25*X.Test+0.25*U.Test))
# error.Test <- 0.25*rnorm(N)+0.25*U.Test
# Y.Test <- 1+0.25*X.Test+(0.5+0.25*X.Test)*A.Test+error


################################################################################
# Cross-Validation
################################################################################

CrossVar <- function(PARAMETER,
                     ss,
                     subset,
                     response,
                     target,
                     perturb,
                     diagonal ){
  
  P.P <- rep(as.numeric(PARAMETER[1]),2)  # put dimension of perturbing variable (which is Z,L)
  P.T <- rep(as.numeric(PARAMETER[2]),2)  # put dimension of target variable (which is A,L)
  Lambda  <- as.numeric(PARAMETER[3])     # regularization parameter
  
  RISK <- rep(0,NumCV)
  
  for(cv in 1:NumCV){
    
    CV.Split.Index <- list()
    CV.Split.Index[[1]] <- intersect( SS.Index[[ss]][ -CV.Index[[cv]] ],subset)
    CV.Split.Index[[2]] <- intersect( SS.Index[[ss]][ CV.Index[[cv]] ] ,subset)
    
    response.CV  <- list()
    target.CV    <- list()
    perturb.CV   <- list()
    diagonal.CV   <- list()
    
    for(cvest in 1:2){
      response.CV[[cvest]] <- response[CV.Split.Index[[cvest]]]
      target.CV[[cvest]]   <- target[CV.Split.Index[[cvest]],]
      perturb.CV[[cvest]]  <- perturb[CV.Split.Index[[cvest]],]
      diagonal.CV[[cvest]] <- diagonal[CV.Split.Index[[cvest]]]
    }
    
    CV.result <- 
      FT_PMMR_CV(Y.Train       =response.CV[[1]],
                 Perturb.Train =perturb.CV[[1]],
                 Target.Train  =target.CV[[1]],
                 Diagonal.Train=diagonal.CV[[1]],
                 Y.Valid       =response.CV[[2]],
                 Perturb.Valid =perturb.CV[[2]],
                 Target.Valid  =target.CV[[2]],
                 Diagonal.Valid=diagonal.CV[[2]],
                 Perturb.bw    =exp(P.P),
                 Target.bw     =exp(P.T),
                 lambda        =exp(Lambda),
                 NV            =FALSE)
    
    RISK[c(cv)] <- c(CV.result$Vstat)
  }
  
  return( mean(RISK) )
}



################################################################################
# Cross-fitting and Cross Validation Sets
################################################################################


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

X.pos <- which( substr(colnames(Data),1,1)=="X" )
Z.pos <- which( substr(colnames(Data),1,1)=="Z" )
A.pos <- which( substr(colnames(Data),1,1)=="A" )
Y.pos <- which( substr(colnames(Data),1,1)=="Y" )

AX <- cbind(A,X)
ZX <- cbind(Z,X)

Y.MM  <- list()
A.MM  <- list()
X.MM  <- list()
Z.MM  <- list()
AX.MM <- list()
ZX.MM <- list()


for(ss in 1:CF){
  
  Y.MM [[ss]]  <- Y[SS.Index[[ss]] ]
  A.MM [[ss]]  <- A[SS.Index[[ss]] ]
  X.MM [[ss]]  <- X[SS.Index[[ss]] ]
  Z.MM [[ss]]  <- Z[SS.Index[[ss]] ]
  AX.MM[[ss]]  <- cbind(A.MM[[ss]],X.MM[[ss]])
  ZX.MM[[ss]]  <- cbind(Z.MM[[ss]],X.MM[[ss]])
  
}

# ###################################################
# # Find opt hyperpara for PMMR: Commented out because too slow
# ###################################################
# 
# Opt.Para.h <- list()
# for(ss in 1:CF){
#   
#   CV.Curve <- apply(Para.Grid,
#                     1,
#                     FUN=function(vv){ CrossVar(vv,
#                                                ss,
#                                                subset=1:N,
#                                                response=Y,
#                                                target=AX,
#                                                perturb=ZX,
#                                                diagonal=rep(1,N)) })
#   
#   Opt.Para.h[[ss]] <-
#     as.numeric(Para.Grid[which.min(CV.Curve),])
# }

# Comment these two lines if the hyperparameters are chosen from cross-validation
Opt.Para.h <- list()
Opt.Para.h[[1]] <- Opt.Para.h[[2]] <- c(0,1,-2)

h.MM <- list()
h.predict <- list()

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

################################################################################
# Summary
################################################################################

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

ATT.Estimate <- h.A1.hat.CF-h.A0.hat.CF

plot(ATT.true,ATT.Estimate); abline(a=0,b=1,col=2,lwd=2)

mean(ATT.true)
mean(ATT.Estimate)










