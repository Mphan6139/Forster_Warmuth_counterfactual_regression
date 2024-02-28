library(MASS)
library(splines)
library(latex2exp)
source('toSource.R')


expit <- function(x){ exp(x)/(1+exp(x)) };

logit <- function(x){ log(x/(1-x)) }

n <- 4*2000; nsim <- 500;

rateseq <- seq(0.1,0.5,by=0.05);



set.seed(1)
for (D in c(4)){
  res <- cate_hat <- res2 <- NULL
  for (rate in rateseq){
    
    print(rate)
    
    cate_hat_mat <- data.frame(matrix(nrow=nsim,ncol=6))
    res <- data.frame(matrix(nrow=nsim,ncol=6))
    colnames(res) <- c("forster_poly", "ls_poly","forster_ns","ls_ns", "forster_bs", "ls_bs")
    for (i in 1:nsim){
      if(i %%100==0){print(i)}
      s <- sort(rep(1:4,n/4));
      x <-  runif(n, -1, 1);
      ps <- 0.1 + 0.8*(x>0)
      mu0 <- (x <= -.5)*0.5*(x+2)^2 + (x/2+0.875)*(x>-1/2 & x<0) +
        (x>0 & x<.5)*(-5*(x-0.2)^2 +1.075) + (x>.5)*(x+0.125);
      mu1 <- mu0;
      tau <- 0
      a <- rbinom(n,1,ps);
      y <- a*mu1 + (1-a)*mu0 + rnorm(n,sd=(.2-.1*cos(2*pi*x)))
      
      ## estimate nuisance functions
      
      pihat <- expit( logit(ps) + rnorm(n,mean=1/(n/4)^rate,sd=1/(n/4)^rate))
      
      mu1hat <- predict(smooth.spline(x[a==1 & s==2],y[a==1 & s==2], cv = T),x)$y
      mu0hat <- predict(smooth.spline(x[a==0 & s==2],y[a==0 & s==2], cv = T),x)$y
      
      ## construct estimators
      
      plugin <- mu1hat-mu0hat
      
      pseudo <- ((a-pihat)/(pihat*(1-pihat)))*(y-a*mu1hat-(1-a)*mu0hat) + mu1hat-mu0hat
      
      
      
      tau_forster_poly = series_df(x[s==3], pseudo[s==3], x[s==4],df=D, type = "forster", basis_type = "poly")
      tau_forster_ns = series_df(x[s==3], pseudo[s==3], x[s==4],df=D, type = "forster", basis_type = "ns")
      tau_forster_bs = series_df(x[s==3], pseudo[s==3], x[s==4],df=D, type = "forster", basis_type = "bs")
      
      tau_ls_poly = series_df(x[s==3], pseudo[s==3], x[s==4], df=D, type = "ls", basis_type = "poly")
      tau_ls_ns = series_df(x[s==3], pseudo[s==3], x[s==4], df=D, type = "ls", basis_type = "ns")
      tau_ls_bs = series_df(x[s==3], pseudo[s==3], x[s==4], df=D, type = "ls", basis_type = "bs")
      
      
      res$forster_poly[i] <- (n/4)*mean((tau_forster_poly[[1]]-tau)^2)
      res$forster_ns[i] <- (n/4)*mean((tau_forster_ns[[1]]-tau)^2)
      res$forster_bs[i] <- (n/4)*mean((tau_forster_bs[[1]]-tau)^2)
      res$ls_poly[i] <-(n/4)*mean((tau_ls_poly[[1]]-tau)^2)
      res$ls_ns[i] <- (n/4)*mean((tau_ls_ns[[1]] -tau)^2)
      res$ls_bs[i] <- (n/4)*mean((tau_ls_bs[[1]] -tau)^2)
      
      cate_hat_mat$forster_poly[i] <- mean(tau_forster_poly[[1]] )
      cate_hat_mat$forster_ns[i] <- mean(tau_forster_ns[[1]] )
      cate_hat_mat$ls_poly[i] <-mean(tau_ls_poly[[1]] )
      cate_hat_mat$ls_ns[i] <- mean(tau_ls_ns[[1]] )
    }
    
    
    res2 <- rbind(res2, c(rate, apply(res,2,mean)))
    cate_hat<- rbind(cate_hat, c(rate, apply(cate_hat_mat,2,mean)))
  }
}  

save(res2, file = paste0('uniform_forster_ls_',D,'_',n/4,'.RData'))


nn <- ncol(res2)
pdf(file=paste0("uniform_forster_ls_df_",D,"n=",n/4,".pdf"))
layout(matrix(c(1,2),nrow=1), width=c(4,1.5)) 
par(mar=c(5,4,4,0)) #No margin on the right side
matplot(rateseq, res2[,-1], type = "o", ylab = "n* MSE", 
        xlab = TeX(r'($\hat{\pi}$ convergence rate ($\alpha$ in RMSE $n^{-\alpha})$)'))
title(main=paste0("Forster vs ls, df=",D,"n=",n/4))
par(mar=c(5,0,4,2)) #No margin on the left side
plot(c(0,1),type="n", axes=F, xlab="", ylab="")

ltyp <- rep(1:5, times=nn/5, each=1)
cols <- rep(1:6, times=nn/6, each=1)
legend("top", pch = as.character(1:nn), colnames(res2)[-1], lty = ltyp, col=cols,cex=0.8)
dev.off()



pdf(file=paste0("uniform_forster_ls_df_",D,"n=",n/4,"_splines.pdf"))
layout(matrix(c(1,2),nrow=1), width=c(4,1.5)) 
par(mar=c(5,4,4,0)) #No margin on the right side
matplot(rateseq, res2[,4:7], type = "o", ylab = "n* MSE", 
        xlab = TeX(r'($\hat{\pi}$ convergence rate ($\alpha$ in RMSE $n^{-\alpha})$)'))
title(main=paste0("Forster vs ls, df=",D,"n=",n/4, ', with splines'))
par(mar=c(5,0,4,2)) #No margin on the left side
plot(c(0,1),type="n", axes=F, xlab="", ylab="")

ltyp <- rep(1:5, times=nn/5, each=1)
cols <- rep(1:6, times=nn/6, each=1)
legend("top", pch = as.character(1:nn), colnames(res2)[4:7], lty = ltyp, col=cols,cex=0.8)
dev.off()












