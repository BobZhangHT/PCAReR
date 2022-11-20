# this script generates the data

# Covariate generation ----------------------------------------------------
x_gen <- function(d, rho, n=100) {
  library(MASS)
  # d: dimension of covariates
  # rho: covariance
  # n: sample size
  mu = rep(0,d)
  sigmat = matrix(rho,d,d)
  diag(sigmat) = 1
  X = mvrnorm(n,mu,sigmat)
  # X = scale(X,scale = F) # center the data
  return(X)
}

# Response generation -----------------------------------------------------
y_gen = function(X,w,beta,tau=1,std_err=1,type='linear'){
  # X: design matrix
  # w: allocation vector
  # beta: coef
  # tau: treatment effect
  # std_err: residual std error
  n = nrow(X)
  if(type=='linear'){
    gfun = X%*%beta
  }else if(type=='nonlinear1'){
    gfun = exp(X)%*%beta
  }else if(type=='nonlinear2'){
    gfun = exp(2*X%*%beta/d)
  }else if(type=='nonlinear3'){
    gfun = sin(pi*X%*%beta/d)
  }else if(type=='nonlinear4'){
    gfun = (X%*%beta/d)^2
  }
  y = gfun+tau*w+std_err*rnorm(n)
  return(y)
}


