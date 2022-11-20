rm(list=ls())

library(glmnet)
library(fastDummies)
library(latex2exp)
library(ggplot2)

# load functions
src_dir = './src/'
files = list.files(src_dir)
for(file in files){
  source(paste(src_dir,file,sep=''))
}

# load the data from Hill (2011)
load('./data/ucgs_a_10712097_sm0001/data/example.data')
load('./data/ucgs_a_10712097_sm0001/data/sim.data')

# combine imp1 with the treatment effect
imp1 = data.frame(iqsb.36=ihdp$iqsb.36,imp1)
covs.cont.n=c("bw","b.head","preterm","birth.o","nnhealth","momage")
covs.cat.n=c("sex","twin","b.marr","mom.lths","mom.hs",	"mom.scoll","cig","first","booze","drugs","work.dur","prenatal","ark","ein","har","mia","pen","tex","was")
p=length(c(covs.cont.n,covs.cat.n))
imp1 = imp1[,c('iqsb.36','treat',covs.cont.n,covs.cat.n)]

for(cov_name in covs.cat.n){
  cat(cov_name,unique(imp1[,cov_name]),'\n')
}
imp1[,'first'] = ifelse(imp1[,'first']==2,1,0)
imp1=na.omit(imp1) # removing the missing values
ihdp = imp1

# Data Preprocessing ------------------------------------------------------

# this data set can be found on the repo
# https://github.com/AMLab-Amsterdam/CEVAE

colnames(ihdp) = c('y','treatment',
                   paste('x',1:(ncol(ihdp)-2),sep=''))
# standardize the continuos covariates
#ihdp[,paste('x',1:6,sep='')] = scale(ihdp[,paste('x',1:6,sep='')])

# standardize all comvariates
ihdp[,1:27] = scale(ihdp[,1:27])

# Simulation --------------------------------------------------------------
# fit a linear model to simulate the real process
# all variables
set.seed(2020)
vars = paste(paste('x',1:(ncol(ihdp)-2),sep=''),collapse = '+')
# formula
fml = as.formula(paste('~ (',vars,') ^2 - 1'))
design_mat = data.frame(treatment=ihdp$treatment,
                        model.matrix(fml,data=ihdp))
design_mat = scale(design_mat)
cv_output = cv.glmnet(x=as.matrix(design_mat),y=ihdp$y,
                      family = 'gaussian')
best_lam <- cv_output$lambda.min
best_lam
lm_fit = glmnet(x=as.matrix(design_mat),y=ihdp$y,
                lambda=best_lam,
                family = 'gaussian')
tau = coef(lm_fit)[2]
sum(coef(lm_fit)!=0)

# a function to generate the response
ihdp_y_gen = function(X,w){
  XX = cbind(w,X)
  yhat = predict(lm_fit,XX)
  yhat
}

# We balance all samples
n = nrow(ihdp)
X = as.matrix(model.matrix(fml,data=ihdp))[1:n,]
X = scale(X)

save_nm = '0608-ihdp-all2ways'

# Determine K
X_svd = svd(X)
scree_df = data.frame(
  index = 1:dim(X)[2],
  var_explained=X_svd$d^2/sum(X_svd$d^2)
)

cumsum_vars = cumsum(X_svd$d^2)
cumsum_vars = cumsum_vars/cumsum_vars[length(cumsum_vars)]

scree_plt = ggplot(scree_df,aes(x=index,y=var_explained, group=1))+
  geom_line()+
  xlab('Index')+
  ylab('Variance Explained')+
  geom_vline(aes(xintercept=sum(X_svd$d^2>mean(X_svd$d^2)),
                linetype='solid'),
             color='darkred',
             show.legend = T,#'Kaiser',
             key_glyph = "path",
             lwd=0.5)+
  geom_vline(aes(xintercept=sum(cumsum_vars<0.5)+1,
                linetype='dashed'),
             color='darkblue',
             show.legend = T,
             key_glyph = "path",
             lwd=0.5)+
  scale_linetype_manual("", 
                      values = c('solid','dashed'),
                      guide = guide_legend(override.aes = list(color = c("darkblue", 
                                                                         "darkred"))),
                      labels = unname(TeX(c('$\\gamma_k=0.5$',
                                               'Kaiser')))) +
  #scale_linetype_discrete(values=c('solid','dashed'))+
  theme_bw()+
  theme(legend.position="top")


# given X, perform all randomization methods for 1000 times
reps = 1000
n_budget = 10000
pa = 0.05
var_ratio = 0.5

n0 = ceiling(n/2)
n1 = n-n0

# parallel computation
library(abind)
library(foreach)
library(doRNG)
library(doSNOW)

cores = 50

# create the folder
save_dir_folder = file.path('./save',save_nm)
if(!file.exists(save_dir_folder)){
  dir.create(save_dir_folder,recursive = T)
}

# generate allocation matrix with the dimension of n x reps
save_dir_rand = file.path('./save',save_nm,'cr.rdata')
if(!file.exists(save_dir_rand)){
  set.seed(2020)
  
  start_time = Sys.time()
  W_rand = c()
  for(i in 1:reps){
    W_rand = cbind(W_rand,sample(rep(c(0,1),c(n0,n1)),n,replace = F))
  }
  end_time = Sys.time()
  time_rand = as.numeric(difftime(end_time, start_time, units="secs"))
  
  save(W_rand,time_rand,file = save_dir_rand)
}else{
  load(file = save_dir_rand)
}

save_dir_rer = file.path('./save',save_nm,'rer.rdata')
if(!file.exists(save_dir_rer)){
  set.seed(2020)
  ii_vec_rer = c()
  W_rer = c()
  
  start_time = Sys.time()
  n = nrow(X)
  d = ncol(X)
  a = qchisq(pa,d)
  
  n0 = ceiling(n/2)
  n1 = n-n0
  
  for(i in 1:reps){
    
    w = sample(rep(c(0,1),c(n0,n1)),n,F) # generate allocation
    mdist = maha_dist(X,w)
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
    ii_vec_rer = c(ii_vec_rer,ii)
    
    indicator = ii>=n_budget
    if(indicator == T){
      W_rer = cbind(W_rer,as.numeric(best_w))
    }else{
      W_rer = cbind(W_rer,as.numeric(w))
    }
  }
  end_time = Sys.time()
  time_rer = as.numeric(difftime(end_time, start_time, units="secs"))
  
  save(W_rer,time_rer,
       ii_vec_rer,
       file = save_dir_rer)
}else{
  load(file = save_dir_rer)
}

save_dir_pcarer = file.path('./save',save_nm,'pcarer.rdata')
if(!file.exists(save_dir_pcarer)){
  ii_vec_pcarer = c()
  set.seed(2020)
  
  W_pcarer = c()
  
  start_time = Sys.time()
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
  
  for(i in 1:reps){
    
    w = sample(rep(c(0,1),c(n0,n1)),n,F) # generate allocation
    mdist = maha_dist(Zk,w)
    ii = 1
    best_w = w
    best_mdist = mdist
    while (mdist>a_pca) {
      #w = sample(rep(c(0,1),c(n/2,n/2)),n,F)
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
    ii_vec_pcarer = c(ii_vec_pcarer,ii)
    indicator = ii>=n_budget
    if(indicator == T){
      W_pcarer = cbind(W_pcarer,as.numeric(best_w))
    }else{
      W_pcarer = cbind(W_pcarer,as.numeric(w))
    }
  }
  end_time = Sys.time()
  time_pcarer = as.numeric(difftime(end_time, start_time, units="secs"))

  save(W_pcarer,time_pcarer,
       ii_vec_pcarer,
       file = save_dir_pcarer)
}else{
  load(file = save_dir_pcarer)
}

save_dir_ridgerer = file.path('./save',save_nm,'ridgerer.rdata')
if(!file.exists(save_dir_ridgerer)){
  set.seed(2020)
  ii_vec_ridgerer = c()
  W_ridgerer = c()
  
  start_time = Sys.time()
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
  
  for(i in 1:reps){
    
    w = sample(rep(c(0,1),c(n0,n1)),n,F) 
    mdist = maha_dist(X,w,lambda)
    
    ii = 1
    best_w = w
    best_mdist = mdist
    while (mdist>a) {
      #w = sample(rep(c(0,1),c(n/2,n/2)),n,F)
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
    ii_vec_ridgerer = c(ii_vec_ridgerer,ii)
    indicator = ii>=n_budget
    if(indicator == T){
      W_ridgerer = cbind(W_ridgerer,as.numeric(best_w))
    }else{
      W_ridgerer = cbind(W_ridgerer,as.numeric(w))
    }
  }
  end_time = Sys.time()
  time_ridgerer = as.numeric(difftime(end_time, start_time, units="secs"))
  
  save(W_ridgerer,time_ridgerer,
       ii_vec_ridgerer,
       file = save_dir_ridgerer)
}else{
  load(file = save_dir_ridgerer)
}

# save.image('./save/ihdp_workspace.rdata')

# Covariate Analysis ------------------------------------------------------
library(lattice)
library(ggplot2)
library(reshape2)
library(scales)

# covariate var reduction
cov_diff_rand = apply(W_rand, MARGIN=2, FUN=mean_diff, X=X)
cov_diff_rer = apply(W_rer, MARGIN=2, FUN=mean_diff, X=X)
cov_diff_ridgerer = apply(W_ridgerer, MARGIN=2, FUN=mean_diff, X=X)
cov_diff_pcarer = apply(W_pcarer, MARGIN=2, FUN=mean_diff, X=X)

var_mdiff_rand = apply(cov_diff_rand, 2, var)
var_mdiff_rer = apply(cov_diff_rer, 2, var)
var_mdiff_ridgerer = apply(cov_diff_ridgerer, 2, var)
var_mdiff_pcarer = apply(cov_diff_pcarer, 2, var)

r_mdiff_rer = (var_mdiff_rand - var_mdiff_rer)/(var_mdiff_rand+1e-30)
r_mdiff_ridgerer = (var_mdiff_rand - var_mdiff_ridgerer)/(var_mdiff_rand+1e-30)
r_mdiff_pcarer = (var_mdiff_rand - var_mdiff_pcarer)/(var_mdiff_rand+1e-30)
r_mdiff_mat = as.matrix(data.frame(ReR=r_mdiff_rer,
                                   RidgeReR=r_mdiff_ridgerer,
                                   PCAReR=r_mdiff_pcarer))

r_mdiff_all_mat = data.frame(Methods = c('ReR','RidgeReR','PCAReR'),
                             value=colMeans(r_mdiff_mat))
r_mdiff_all_mat$Methods = factor(r_mdiff_all_mat$Methods,
                                 levels = c('ReR','RidgeReR','PCAReR'))
rownames(r_mdiff_all_mat) = NULL

r_mdiff_mat_plt = r_mdiff_mat
rownames(r_mdiff_mat_plt) = NULL
colnames(r_mdiff_mat_plt) = c('ReR','Ridge-ReR','PCA-ReR')

nmlist <- list(
  layout.heights = list(
    top.padding = 0,
    main.key.padding = 0,
    key.axis.padding = 0,
    axis.xlab.padding = 0,
    xlab.key.padding = 0,
    key.sub.padding = 0,
    bottom.padding = 0
  ),
  layout.widths = list(
    left.padding = 0,
    key.ylab.padding = 0,
    ylab.axis.padding = 0,
    axis.key.padding = 0,
    right.padding = 0.6
  ),
  axis.components = list(
    top = list(pad1 = 0, pad2 = 0), # padding above top axis
    right = list(pad1 = .2, pad2 = .2)
  )
)

plt_mdiff = levelplot(r_mdiff_mat_plt[1:25,c(3,2,1)],
                   xlab=list('Covariate Index',cex=.75),
                   ylab='',
                   aspect = 0.3,
                   #xlab=TeX('$X$'),ylab='',#TeX('$d$'),
                   #at=seq(0,0.5,.1),
                   #main=TeX("$r_{\}$"))
                   par.settings=nmlist,
                   main='')

plt_mdiff

pdf(file.path('./save',save_nm,'ihdp_cov.pdf'),
    height=2,width=6)
plt_mdiff
dev.off()

# $\tau$ Analysis  --------------------------------------------------------

# analyze tau 
tau_rand = apply(W_rand, MARGIN = 2, FUN = function(w){
  y = ihdp_y_gen(X,w)
  tau_est(y, w)
})

tau_rer = apply(W_rer, MARGIN = 2, FUN = function(w){
  y = ihdp_y_gen(X,w)
  tau_est(y, w)
})

tau_pcarer = apply(W_pcarer, MARGIN = 2, FUN = function(w){
  y = ihdp_y_gen(X,w)
  tau_est(y, w)
})

tau_ridgerer = apply(W_ridgerer, MARGIN = 2, FUN = function(w){
  y = ihdp_y_gen(X,w)
  tau_est(y, w)
})

tau_df = data.frame(c(tau_rer,tau_ridgerer,tau_pcarer))
colnames(tau_df) = 'value'
tau_df['Methods'] = factor(rep(c('ReR','Ridge-ReR','PCA-ReR'),rep(reps,3)),
                          levels=c('ReR','Ridge-ReR','PCA-ReR'))
library(latex2exp)

denplt<-ggplot(tau_df, aes(x=value, 
                           color=Methods,
                           fill=Methods,
                           linetype=Methods
                           )) +
  geom_density(alpha=0.4)+
  geom_vline(aes(xintercept=tau),
             linetype="dashed")+
  xlab(TeX('$\\hat{\\tau}$'))+
  ylab('Density') + 
  #scale_color_grey(start = 0.2,
  #                 end = 0.6) + 
  theme_bw()
ggsave(file.path('./save',save_nm,'ihdp_denplot.pdf'),
       plot=denplt,width = 5,height = 4)
denplt
cat('Density Plot Saved!\n')

# save as a table
tau_df = data.frame(CR = tau_rand,
                    ReR = tau_rer,
                    RidgeReR = tau_ridgerer,
                    PCAReR = tau_pcarer)
tau_df = apply(tau_df,2,function(x){
  mean((x-tau)^2)
})
tau_df = 1- tau_df[2:4]/tau_df[1]
tau_df = data.frame(Methods = c('ReR','RidgeReR','PCAReR'),
                    value=tau_df)
tau_df$Methods = factor(tau_df$Methods,
                        levels = c('ReR','RidgeReR','PCAReR'))
rownames(tau_df) = NULL

df_time_mdiff_tau = data.frame(Methods=c('ReR','RidgeReR','PCAReR'),
                               mdiff=r_mdiff_all_mat[,2],
                               tau=tau_df[,2],
                               time=c(mean(time_rer),
                                      mean(time_ridgerer),
                                      mean(time_pcarer)))
df_time_mdiff_tau
write.csv(df_time_mdiff_tau,
          file.path('./save',save_nm,'ihdp_table.csv'),
          row.names = F)

ggsave(file.path('./save',save_nm,'ihdp_scree_plot.pdf'),
       plot=scree_plt,width = 5,height = 4)
scree_plt
cat('Scree Plot Saved!\n')
