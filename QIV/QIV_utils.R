###
single_gen_QIV = function(n){
  ###
  # Simulation Parameters 10/10
  # p_z = 0.5
  # pi_1 = 0.75
  # pi_0 = 0.25
  # phi_1 = 1
  # psi = 4
  ###
  
  # Covariates:
  #x1 <- runif(n,0,1)
  #x2 <- rnorm(n,0,2)
  x <- rnorm(n,0,1)
  u <- rnorm(n,0,1)
  
  ### Correlated
  # ux <- MASS::mvrnorm(n, mu = c(0,0), Sigma = matrix(c(1,0.5,0.5,1),2,2))
  ###
  
  # IV:
  # Restriction: Independence from U
  p_z = 0.5
  z <- rbinom(n,1,prob = p_z)
  
  # Treatment:
  ps <- function(z,u){
    p_0 <- 0.25
    alpha_z <- log(3)
    alpha_u <-function(u){0.1*(u>0)}
    return(p_0*exp(alpha_z*z + alpha_u(u)))
  }
  # generate p_0 as a function of x (expit) for potentially nonlinear x
  ps_ = ps(z,u)
  a <- rbinom(n,1,ps_);
  
  
  # Outcome
  # Restrictions: Exclusion restriction and no A,U interaction.
  p_2 <- function(Z,U,A,X){
    #beta_a <- function(u){1*(u>0)}
    beta_a <- function(u){4}
    
    beta_u <- function(u){0.5*u}
    beta_z <- 0.5
    beta_x <- 0.25
    
    return(beta_a(U)*A + beta_u(U) + beta_z*Z + beta_x*X)
  }
  error <- rnorm(n,sd=0.2)
  y <- p_2(Z=z,U=u,A=a,X=x) + error
  
  test_6_1 <- mean(p_2(Z=1,U=u,A=1,X=x)-p_2(Z=1,U=u,A=0,X=x)) ==
    mean(p_2(Z=0,U=u,A=1,X=x)-p_2(Z=0,U=u,A=0,X=x))
  
  u_s <- rnorm(n,0,1)
  test_6_2 <- mean(p_2(Z=1,U=u_s,A=1,X=x)-p_2(Z=1,U=u_s,A=0,X=x)) ==
    mean(p_2(Z=1,U=u,  A=1,X=x)-p_2(Z=1,U=u,  A=0,X=x))
  test_6_3 <- mean(p_2(Z=0,U=u_s,A=1,X=x)-p_2(Z=0,U=u_s,A=0,X=x)) ==
    mean(p_2(Z=0,U=u,  A=1,X=x)-p_2(Z=0,U=u,  A=0,X=x))
  
  if(!(test_6_1&test_6_2&test_6_3)){
    stop("assumptions violated!")
  }
  
  
  PHI_1 = p_2(Z=1,U=u,A=1,X=x) - p_2(Z=0,U=u,A=1,X=x)
  delta_a = ps(z=1,u) -ps(z=0,u)
  delta_y = p_2(Z=1,U=u,A=a,X=x) - p_2(Z=0,U=u,A=a,X=x)
  PSI = (delta_y-PHI_1)/delta_a
  mu1z <- p_2(Z=z,U=u,A=rep(1,n),X=x) + error
  mu0z <- p_2(Z=z,U=u,A=rep(0,n),X=x) + error
  ATT <-  sum(a*mu1z-a*mu0z)/sum(a) 
  df <- data.frame("X"=x,"Y"=y,"Z"=z,"A"=a)
  
  return(list(data = df, att = ATT,phi_1 = PHI_1,psi = PSI))
}
