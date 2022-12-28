# update 0220: accommondate unequal sample size

maha_dist = function(X,w,lambda=0){
  # note that the input data matrix should be centered
  library(MASS)
  n1 = sum(w)
  n0 = sum(1-w)
  n = n0+n1
  xbar1 = t(X)%*%w/n1
  xbar0 = t(X)%*%(1-w)/n0
  delta = xbar1-xbar0
  # cov_mat = 4/(n*(n-1))*t(X)%*%X
  cov_mat = (n/(n1*n0*(n-1)))*t(X)%*%X
  # use ginv other than solve to prevent some numerical issues
  mdist = as.numeric(t(delta)%*%ginv(cov_mat+lambda*diag(nrow(cov_mat)))%*%delta)
  return(mdist)
}

# rerandomization
ReR <- function(pa, X, 
                n_budget=1000,
                seed=2020) {
  # pa: threshold probability
  # X: input data matrix
  set.seed(seed)
  
  n = nrow(X)
  d = ncol(X)
  a = qchisq(pa,d)
  
  n0 = ceiling(n/2)
  n1 = n-n0
  w = sample(rep(c(0,1),c(n0,n1)),n,F) # generate allocation
  mdist = ifelse(n<d,n-1,maha_dist(X,w))
  
  # introduce n_budget to avoid the extreme case 
  # when a is too restrictive for rerandomization
  
  ii = 1
  best_w = w
  best_mdist = mdist
  while (mdist>a) {
    w = sample(rep(c(0,1),c(n0,n1)),n,F)
    mdist = maha_dist(X,w)
    if(best_mdist>mdist){
      best_w = w
      best_mdist = mdist
    }
    ii = ii + 1
    if(ii>=n_budget){
      break
    }
  }
  
  indicator = ii>=n_budget
  if(indicator == T){
    return_w = best_w
  }else{
    return_w = w
  }

  return(list(w=return_w,
              ii=ii,
              indicator=indicator,
              a=a))
}

# PCA rerandomization
PCA_ReR = function(pa,X,var_ratio=0.95, 
                   n_budget=1000,seed=2020){
  # pa: threshold prob
  # X: centered design matrix
  # var_ratio: variance ratio to select k
  #            if var_ratio='Kasier', 
  #            we use the Kasier rule
  
  set.seed(seed)
  
  n = nrow(X)
  d = ncol(X)
  
  X_svd = svd(X)
  if (var_ratio=='Kaiser'){
    k = sum(X_svd$d^2>mean(X_svd$d^2))
  }else{
    cumsum_vars = cumsum(X_svd$d^2)
    cumsum_vars = cumsum_vars/cumsum_vars[length(cumsum_vars)]
    k = sum(cumsum_vars<var_ratio)+1
  }
  
  if(k==1){
    Zk = X_svd$u[,1]*X_svd$d[1]
  }else{
    Zk = X_svd$u[,1:k]%*%diag(X_svd$d[1:k])
  }
  
  a_pca = qchisq(pa,k)
  
  n0 = ceiling(n/2)
  n1 = n-n0
  w = sample(rep(c(0,1),c(n0,n1)),n,F) # generate allocation
  mdist = maha_dist(Zk,w)
  
  ii = 1
  best_w = w
  best_mdist = mdist
  while (mdist>a_pca) {
    w = sample(rep(c(0,1),c(n0,n1)),n,F) 
    mdist = maha_dist(Zk,w)
    if(best_mdist>mdist){
      best_w = w
      best_mdist = mdist
    }
    ii = ii + 1
    if(ii>=n_budget){
      break
    }
  }
  
  indicator = ii>=n_budget
  if(indicator == T){
    return_w = best_w
  }else{
    return_w = w
  }
  
  return(list(k=k,w=return_w,
              indicator=indicator,
              ii=ii,
              a=a_pca))
}

# Ridge rerandomization
# the function to find the quantile of the ridge regression estimator
find_q = function(lambda,lambdas,pa=0.05,ub=NULL){
  d = length(lambdas)
  f = function(t,q){
    tmp1 = sin(0.5*(-t*q+sum(atan(lambdas/(lambda+lambdas+1e-8)*t))))
    tmp2 = t*prod(1+(lambdas/(lambda+lambdas+1e-8))^2*t^2)^(1/4)
    tmp1/tmp2
  }
  vf = Vectorize(f,vectorize.args = 't')
  Fq = function(q){
    tryCatch({
      0.5-(1/pi)*integrate(vf,lower = 0,
                           upper = Inf,q=q,
                           subdivisions=5000)$value
    },error=function(e){
      xi = 1e-4
      U = (xi*pi*(d/2)*prod(lambdas/(lambda+lambdas+1e-8)))^(2/d)
      0.5-(1/pi)*integrate(vf,lower = 0,
                           upper = U,q=q,
                           subdivisions=5000)$value
    })
  }
  fn = function(q,pa){
    (Fq(q)-pa)^2
  }
  # fix a bug when seq_ub < 0.2
  seq_ub = ifelse(is.null(ub),d,ub)
  if(seq_ub>0.2){
    init_seq = seq(0.1,seq_ub,0.1)
  }else{
    init_seq = seq(0,seq_ub,length=100)
  }
  init_idx = which.min((sapply(init_seq, Fq)-pa)^2)
  result = nlminb(start=init_seq[init_idx], objective=fn,
                  lower=0, pa=pa)$par
  result
}

find_lambda = function(sigmat,lambdas,V,pa=0.05,
                       delta=1e-2,eps=1e-4,
                       n_search=10,
                       nn=1000){
  
  # sigmat: covariance matrix
  # lambdas: the eigenvalues of sigmat
  # V: the eigenvectors of sigmat
  # pa, delta, eps: the parameters of 'procedure fore finding a desirable lambda>0'
  # n_search: the budget of iterations to search lambda
  # nn: monte carlo budget for estimation
  
  d = nrow(V)
  # initialization
  lambda = 0 
  Lambda = c()
  dk_mat = c()
  a_seq = c()
  Z = matrix(rnorm(nn*d),nn,d)
  
  # first calculation
  a = find_q(lambda,lambdas,pa)
  a_plus = find_q(lambda+delta,lambdas,pa,a)
  
  i = 0
  while (abs((lambda+delta)*a_plus-lambda*a)>eps) {
    #cat(abs((lambda+delta)*a_plus-lambda*a))
    #cat('\n')
  
    lambda = lambda + delta
    a = a_plus
    a_plus = find_q(lambda+delta,lambdas,pa,a)
    
    M_lambda = Z^2%*%(lambdas/(lambda+lambdas))
    indicator = M_lambda<=a
    dk = t(indicator)%*%Z^2/(sum(indicator)+eps)
    vk = diag(V%*%diag(lambdas*as.vector(dk))%*%t(V))/(diag(sigmat)+1e-8)
    if(mean(vk) < pchisq(qchisq(pa,df=d),d+2)/pa){
      Lambda = c(Lambda,lambda)
      dk_mat = rbind(dk_mat,dk)
      a_seq = c(a_seq,a)
    }
    
    i = i + 1
    if(i>n_search){
      break
    }
  }
  
  if(length(Lambda)==0){
    lambda = 0
    a = find_q(0,lambdas,pa)
  }else{
    ck = lambdas^2/sum(lambdas^2)
    objs = dk_mat^2%*%ck - (dk_mat%*%ck)^2
    lambda = Lambda[which.min(objs)]
    a = a_seq[which.min(objs)]
  }
  
  return(list(lambda=lambda,a=a))
}

RidgeReR <- function(pa, X,
                     n_budget=1000,
                     seed=2020) {
  
  set.seed(seed)
  
  # pa: threshold probability
  # X: input data matrix
  n = nrow(X)
  d = ncol(X)
  
  sigmat = cov(X)*4/n
  eig_decomp = eigen(sigmat)
  lambdas = eig_decomp$values
  V = eig_decomp$vectors
  
  params = find_lambda(sigmat,lambdas,V,pa)
  a = params$a
  lambda = params$lambda
  
  n0 = ceiling(n/2)
  n1 = n-n0
  w = sample(rep(c(0,1),c(n0,n1)),n,F) 
  mdist = maha_dist(X,w,lambda)
  # introduce n_budget to avoid the extreme case 
  # when no allocation is satisfactory
  ii = 1
  best_w = w
  best_mdist = mdist
  while (mdist>a) {
    w = sample(rep(c(0,1),c(n0,n1)),n,F) 
    mdist = maha_dist(X,w,lambda)
    if(best_mdist>mdist){
      best_w = w
      best_mdist = mdist
    }
    ii = ii + 1
    if(ii>=n_budget){
      break
    }
  }
  
  indicator = ii>=n_budget
  if(indicator == T){
    return_w = best_w
  }else{
    return_w = w
  }
  
  return(list(w=return_w,
              indicator=indicator,
              ii=ii,
              a=a))
}
