rm(list=ls())
library(reshape2)
library(forecast)

# input the directory of your files here
folder = './save/...'

r_dim_df = pca_dim_df
r_dim_df[,4] = pca_dim_df[,4]/as.integer(dim_df[,1]*dim_df[,2])
r_dim_df

# covariate balance -------------------------------------------------------
load(paste0(folder,'summary_tb.rdata'))
cov_df = melt(diff_ratio_df,id.vars = c('n','rnd','rho'),
              variable.name = 'scheme')
for(i in 1:4){
  cov_df[,i] = as.factor(cov_df[,i])
}
cov_df[,5] = 1-cov_df[,5]
cov_df

agg_cov_df = aggregate(value~rnd+rho+scheme,data = cov_df, mean)
agg_cov_df
rfmt_cov_df = dcast(data=agg_cov_df,scheme~rho+rnd)
rfmt_cov_df

agg_cov_k_df = aggregate(pca_dim~rnd+rho,data=r_dim_df,mean)
agg_cov_k_df['scheme'] = 'PCAReR-dim'
rfmt_cov_k_df = dcast(data=agg_cov_k_df,scheme~rho+rnd,value.var = 'pca_dim')
rfmt_cov_k_df

rfmt_cov_all_df = rbind(rfmt_cov_df,rfmt_cov_k_df)
rfmt_cov_all_df

# save
write.csv(rbind(rfmt_cov_df,rfmt_cov_k_df),
          file = paste0(folder,'table_cov.csv'),
          row.names = F)


# time --------------------------------------------------------------------
time_aov_df = melt(time_df,id.vars = c('n','rnd','rho'),
                   variable.name = 'scheme')
for(i in 1:4){
  time_aov_df[,i] = as.factor(time_aov_df[,i])
}

agg_time_df = aggregate(value~rnd+n+scheme,data = time_aov_df, mean)
agg_time_df
rfmt_time_df = dcast(data=agg_time_df,scheme~n+rnd)
rfmt_time_df

agg_time_k_df = aggregate(pca_dim~rnd+n,data=r_dim_df,mean)
agg_time_k_df['scheme'] = 'PCAReR-dim'
rfmt_time_k_df = dcast(data=agg_time_k_df,scheme~n+rnd,value.var = 'pca_dim')
rfmt_time_k_df

write.csv(rbind(rfmt_time_df,rfmt_time_k_df),
          file = paste0(folder,'table_time.csv'),
          row.names = F)

# tau ---------------------------------------------------------------------
tau_df = melt(mse_df,id.vars = colnames(mse_df)[1:6],
              variable.name = 'scheme')
for(i in 1:7){
  tau_df[,i] = as.factor(tau_df[,i])
}
tau_df[,8] = 1-tau_df[,8]

agg_tau_df = aggregate(value~rnd+n+ytype+scheme,data = tau_df, mean)
agg_tau_df
rfmt_tau_df = dcast(data=agg_tau_df,ytype+scheme~n+rnd)
rfmt_tau_df

write.csv(rfmt_tau_df,
          file = paste0(folder,'table_tau.csv'),
          row.names = F)
