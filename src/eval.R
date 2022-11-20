mean_diff = function(X,w){
  n1 = sum(w)
  n0 = sum(1-w)
  xbar1 = t(X)%*%w/n1
  xbar0 = t(X)%*%(1-w)/n0
  delta = xbar1-xbar0
  return(delta)
}

tau_est = function(y,w){
  mean(y[w==1])-mean(y[w==0])
}

infer_fun <- function(M, X, y, method, w, pa, var_ratio, tau) {
  library(ri)
  n = nrow(X)
  perms = replicate(M,{
    if(method == 'Rand'){
      sample(rep(c(0,1),c(n/2,n/2)),n,F)
    }else if(method == 'ReR'){
      ReR(pa,X)
    }else if(method == 'PCA-ReR'){
      PCA_ReR(pa,X,var_ratio)$w
    }else if(method == 'Ridge-ReR'){
      RidgeReR(pa,X)
    }
  })
  
  # confidence interval
  probs = genprob(perms)
  ci_lb = invert.ci(y,w,probs,perms,0.025)
  ci_ub = invert.ci(y,w,probs,perms,0.975)
  
  ci_len = as.numeric(ci_ub - ci_lb)
  ci_cover = as.numeric((tau<=ci_ub)&(tau>=ci_lb))
  
  # p-value
  tau_hat = apply(perms, 2, function(x){
    mean(y[x==1])-mean(y[x==0])
  })
  tau_obs = tau_est(y,w)
  pval = (1+sum(abs(tau_hat)>=abs(tau_obs)))/(1+M)
  
  return(c(ci_len,ci_cover,pval))
}
