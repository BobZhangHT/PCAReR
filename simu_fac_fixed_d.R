rm(list=ls())

# performing the factorial designs simulation following Gutman & Rubins

# load functions
src_dir = './src/'
files = list.files(src_dir)
for(file in files){
  source(paste(src_dir,file,sep=''))
}
# parallel computation
library(foreach)
library(doRNG)
library(doSNOW)
library(abind)

# parameters --------------------------------------------------------------

# experiment
reps = 2000
cores = 50
save_nm = '06-01Simu'

# factors
# known factors
n_seq = c(100,200,500,1000)
d_seq = c(10,50,90,180)
rho_seq = c(0.1,0.5,0.9)

# unknown factors
std_err_seq = c(1,sqrt(0.5))
betatype_seq = c('equal','unequal')
ytype_seq = c('linear','nonlinear1')
# fixed factor
tau = 1

# rerandomization
pa = 0.05
var_ratios = list(0.5,
                  0.7,
                  0.9,
                  'Kaiser')

var_ratios_nm = c('PCAReR-5','PCAReR-7',
                  'PCAReR-9','PCAReR-K')

n_budget = 1000

params = list(reps=reps,cores=cores,
              n_seq=n_seq,d_seq=d_seq,
              rho_seq=rho_seq,
              std_err_seq=std_err_seq,
              betatype_seq=betatype_seq,
              ytype_seq=ytype_seq,
              tau=tau,
              pa=pa,
              var_ratios=var_ratios,
              var_ratios_nm=var_ratios_nm)

# combine function for simulation
cfun = function(list1,list2){
  diff_mat = abind(list1$diff_mat,
                   list2$diff_mat,
                   along = 3)
  tau_mat = abind(list1$tau_mat,
                  list2$tau_mat,
                  along = 3)
  time_vec = abind(list1$time_vec,
                   list2$time_vec,
                   along = 2)
  indicator_vec = abind(list1$indicator_vec,
                        list2$indicator_vec,
                        along = 2)
  iter_vec = abind(list1$iter_vec,
                   list2$iter_vec,
                   along = 2)
  k_vec = abind(list1$k_vec,
                list2$k_vec,
                along = 2)
  return(list(diff_mat=diff_mat,
              tau_mat=tau_mat,
              time_vec=time_vec,
              indicator_vec=indicator_vec,
              iter_vec=iter_vec,
              k_vec=k_vec))
}

cfun_xmat = function(array1,array2){
  array_tmp = abind(array1,array2,along = 3)
  return(array_tmp)
}

# Data Generation ---------------------------------------------------------
save_dir_cache = file.path('./save',save_nm,'cache')
if(!file.exists(save_dir_cache)){
  dir.create(save_dir_cache,recursive = T)
}

d_max = max(d_seq)
n_max = max(n_seq)

set.seed(2020)
for(rho in rho_seq){
  save_dir_cache = file.path('./save',save_nm,'cache',paste0('Xmat_','rho',10*rho,'.rdata'))
  if(!file.exists(save_dir_cache)){
    cat('------------- Gen Xmat for rho =',rho,'-------------\n')
    
    cl <- makeCluster(cores) 
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = reps, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    X_array = foreach(r = 1:reps,
                          .combine = cfun_xmat,
                          .options.snow = opts,
                          .options.RNG = 2020)%dorng%{
                            X = x_gen(d_max,rho,n_max)
                            return(X)
                          }
    close(pb)
    stopCluster(cl)
    
    save(X_array,file = save_dir_cache)
  }
}

# Simulations -------------------------------------------------------------
save_dir = file.path('./save',save_nm)
if(!file.exists(save_dir)){
  dir.create(save_dir,recursive = T)
}

known_facs_comb = expand.grid(n=n_seq,d=d_seq,rho=rho_seq)
unknown_facs_comb = expand.grid(std_err=std_err_seq,
                                betatype=betatype_seq,
                                ytype=ytype_seq)

iter_df = indicator_df = time_df = diff_ratio_df = pca_dim_df = mse_df = c()

for(i in 1:nrow(known_facs_comb)){
  n = known_facs_comb[i,1]
  d = known_facs_comb[i,2]
  rho = known_facs_comb[i,3]
  
  cat('-------------------',
      'n',n,
      'd',d,
      'rho',rho,
      '-------------------\n')
  
  # load X_array
  Xmat_dir_cache = file.path('./save',save_nm,'cache',paste0('Xmat_','rho',10*rho,'.rdata'))
  load(file = Xmat_dir_cache)
  
  # save dir
  save_dir_cache = file.path('./save',save_nm,'cache',paste0('n',n,'d',d,'rho',10*rho))
  if(!file.exists(save_dir_cache)){
    dir.create(save_dir_cache,recursive = T)
  }
  
  # parallel computation
  save_dir_cache_file = paste0(save_dir_cache,'/para_data.rdata')
  if(file.exists(save_dir_cache_file)){
    load(save_dir_cache_file)
  }else{
    cat('run new replications.\n')
    
    cl <- makeCluster(cores) 
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = reps, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    para_comput = foreach(r = 1:reps,
                          .combine = cfun,
                          .options.snow = opts,
                          .options.RNG = 2020)%dorng%{
                            
                            # use the subdata from large x
                            X = X_array[1:n,1:d,r]
                            
                            # standardize the covariates
                            X = scale(X)
                            
                            # randomization
                            start_time = Sys.time()
                            w_rand = sample(rep(c(0,1),c(n/2,n/2)),n,F)
                            end_time = Sys.time()
                            time_rand = as.numeric(difftime(end_time, start_time, units="secs"))
                            
                            # rerandomization
                            start_time = Sys.time()
                            results_rer = ReR(pa,X,n_budget)
                            end_time = Sys.time()
                            w_rer = results_rer$w
                            indicator_rer = results_rer$indicator
                            time_rer = as.numeric(difftime(end_time, start_time, units="secs"))
                            iter_rer = results_rer$ii
                            
                            # ridge rerandomization
                            start_time = Sys.time()
                            results_ridgerer = RidgeReR(pa,X,n_budget)
                            end_time = Sys.time()
                            w_ridgerer = results_ridgerer$w
                            indicator_ridgerer = results_ridgerer$indicator
                            time_ridgerer = as.numeric(difftime(end_time, start_time, units="secs"))
                            iter_ridgerer = results_ridgerer$ii
                            
                            # PCA rerandomization
                            time_pcarer_lst = w_pcarer_lst = list()
                            k_lst = indicator_pcarer_lst = iter_pcarer_lst = list()
                            
                            
                            for(i_var in 1:length(var_ratios)){
                              start_time = Sys.time()
                              results_pcarer = PCA_ReR(pa,X,var_ratios[[i_var]],n_budget)
                              end_time = Sys.time()
                              
                              time_pcarer = as.numeric(difftime(end_time, 
                                                                start_time, 
                                                                units="secs"))
                              w_pcarer = results_pcarer$w
                              k = results_pcarer$k
                              indicator_pcarer = results_pcarer$indicator
                              iter_pcarer = results_pcarer$ii
                              
                              # save results into the list
                              time_pcarer_lst[[i_var]] = time_pcarer
                              w_pcarer_lst[[i_var]] = w_pcarer
                              k_lst[[i_var]] = k
                              indicator_pcarer_lst[[i_var]] = indicator_pcarer
                              iter_pcarer_lst[[i_var]] = iter_pcarer
                            }
                            
                            
                            # covariance evaluation
                            diff_pcarer_lst = list()
                            
                            diff_rand = as.numeric(mean_diff(X,w_rand))
                            diff_rer = as.numeric(mean_diff(X,w_rer))
                            diff_ridgerer = as.numeric(mean_diff(X,w_ridgerer))
                            for(i_var in 1:length(var_ratios)){
                              diff_pcarer_lst[[i_var]] = as.numeric(mean_diff(X,w_pcarer_lst[[i_var]]))
                            }
                            diff_pcarer_mat = matrix(unlist(diff_pcarer_lst),d,length(var_ratios))
                            diff_mat = cbind(diff_rand,diff_rer,
                                             diff_ridgerer,diff_pcarer_mat)
                            colnames(diff_mat) = c('Rand','ReR','RidgeReR',var_ratios_nm)
                            
                            
                            # causal effect estimation
                            tau_rand_vec = c()
                            tau_rer_vec = c()
                            tau_ridgerer_vec = c()
                            tau_pcarer_lst = list()
                            for(ii in 1:nrow(unknown_facs_comb)){
                              
                              std_err = unknown_facs_comb[ii,1]
                              betatype = unknown_facs_comb[ii,2]
                              ytype = unknown_facs_comb[ii,3]
                              
                              if(betatype=='equal'){
                                beta = rep(1,d)
                              }else{
                                beta = c(rep(1,d/2),rep(2,d/2))
                              }
                              
                              y_rand = y_gen(X,w_rand,beta=beta,tau=tau,type=ytype,std_err=std_err)
                              tau_rand_vec = c(tau_rand_vec,tau_est(y_rand,w_rand))
                              
                              y_rer = y_gen(X,w_rer,beta=beta,tau=tau,type=ytype,std_err=std_err)
                              tau_rer_vec = c(tau_rer_vec,tau_est(y_rer,w_rer))
                              
                              y_ridgerer = y_gen(X,w_ridgerer,beta=beta,tau=tau,type=ytype,std_err=std_err)
                              tau_ridgerer_vec = c(tau_ridgerer_vec,tau_est(y_ridgerer,w_ridgerer))
                              
                              
                              tau_pcarer_vec = c()
                              for(i_var in 1:length(var_ratios)){
                                y_pcarer = y_gen(X,w_pcarer_lst[[i_var]],
                                                 beta=beta,tau=tau,type=ytype,std_err=std_err)
                                tau_pcarer_vec = c(tau_pcarer_vec,
                                                   tau_est(y_pcarer,
                                                           w_pcarer_lst[[i_var]]))
                              }
                              tau_pcarer_lst[[ii]] = tau_pcarer_vec 
                              
                            }
                            
                            tau_pcarer_mat = matrix(unlist(tau_pcarer_lst),
                                                    length(var_ratios),
                                                    nrow(unknown_facs_comb))
                            
                            tau_mat = cbind(tau_rand_vec,
                                            tau_rer_vec,
                                            tau_ridgerer_vec,
                                            t(tau_pcarer_mat))
                            
                            colnames(tau_mat) = c('Rand','ReR','RidgeReR',var_ratios_nm)
                            
                            
                            time_vec = c(time_rand,time_rer,time_ridgerer,unlist(time_pcarer_lst))
                            names(time_vec) = c('Rand','ReR','RidgeReR',var_ratios_nm)
                            
                            # exhaust the budget?
                            indicator_vec = c(indicator_rer,indicator_ridgerer,unlist(indicator_pcarer_lst))
                            names(indicator_vec) = c('ReR','RidgeReR',var_ratios_nm)
                            
                            # iterations for each method
                            iter_vec = c(iter_rer,iter_ridgerer,unlist(iter_pcarer_lst))
                            names(iter_vec) = c('ReR','RidgeReR',var_ratios_nm)
                            
                            return(list(diff_mat=diff_mat,
                                        tau_mat=tau_mat,
                                        time_vec=time_vec,
                                        indicator_vec=indicator_vec,
                                        iter_vec=iter_vec,
                                        k_vec=unlist(k_lst)))
                            
                          }
    close(pb)
    stopCluster(cl)
    save(para_comput,file = save_dir_cache_file)
  }
  
  # pvr
  diff_array = para_comput$diff_mat
  diff_vec = colMeans(apply(diff_array, c(1,2), var))
  diff_ratio = diff_vec[2:length(diff_vec)]/diff_vec[1]
  diff_ratio_df_tmp = cbind(known_facs_comb[i,],data.frame(t(diff_ratio)))
  
  # time mat
  time_vec = rowMeans(para_comput$time_vec)[1:length(diff_vec)]
  time_df_tmp = cbind(known_facs_comb[i,],data.frame(t(time_vec)))
  
  # indicator
  indicator_vec = rowMeans(para_comput$indicator_vec)
  indicator_df_tmp = cbind(known_facs_comb[i,],data.frame(t(indicator_vec)))
  
  # iter
  iter_vec = rowMeans(para_comput$iter_vec)
  iter_df_tmp = cbind(known_facs_comb[i,],data.frame(t(iter_vec)))
  
  # k
  pca_dim = rowMeans(para_comput$k_vec)
  names(pca_dim) = var_ratios_nm
  pca_dim_df_tmp = cbind(known_facs_comb[i,],data.frame(t(pca_dim)))

  # tau
  tau_array = para_comput$tau_mat
  mse_mat = apply((tau_array-tau)^2, c(1,2), mean)
  mse_ratio = mse_mat[,2:length(diff_vec)]/mse_mat[,1]
  mse_df_tmp = cbind(known_facs_comb[i,],
                     unknown_facs_comb,
                     mse_ratio,
                     row.names = NULL)
  
  diff_ratio_df = rbind(diff_ratio_df,diff_ratio_df_tmp)
  time_df = rbind(time_df,time_df_tmp)
  indicator_df = rbind(indicator_df,indicator_df_tmp)
  iter_df = rbind(iter_df,iter_df_tmp)
  pca_dim_df = rbind(pca_dim_df,pca_dim_df_tmp)
  mse_df = rbind(mse_df,mse_df_tmp)
}

save(diff_ratio_df,
     time_df,pca_dim_df,
     mse_df,indicator_df,
     iter_df,
     file = file.path(save_dir,'summary_tb.rdata'))

save(params,file = file.path(save_dir,'params.rdata'))

