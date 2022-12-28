rm(list=ls())
library(reshape2)
library(forecast)

# this script generates the ANOVA table and other summary tables

folder = './save/06-07Simu/'

load(paste0(folder,'summary_tb.rdata'))
time_df_summary = time_df
iter_df_summary = iter_df

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

reps = 2000 # total replications
n_reps = 200 # size of the group

# Covariate & Time & Iteration
# Load Data  ------------------------------------------------------------------
known_facs_comb = expand.grid(n=n_seq,d=d_seq,rho=rho_seq)

xdiff_df = c()
time_df = c()
iter_df = c()

for(i in 1:nrow(known_facs_comb)){
  n = known_facs_comb[i,1]
  d = known_facs_comb[i,2]
  rho = known_facs_comb[i,3]
  
  cat('-------------------',
      'n',n,
      'd',d,
      'rho',rho,
      '-------------------\n')
  
  save_dir_cache = file.path(folder,'cache',paste0('n',n,'d',d,'rho',10*rho))
  
  # parallel computation
  save_dir_cache_file = paste0(save_dir_cache,'/para_data.rdata')
  load(save_dir_cache_file)
  
  # covariate and time
  xdiff_array = para_comput$diff_mat
  time_array = para_comput$time_vec
  iter_array = para_comput$iter_vec
  
  mm_df_tmp = c()
  mm_df_tmp_t = c()
  mm_df_tmp_it = c()
  for(i_rep in 1:floor(reps/n_reps)){
    # covariate
    sub_mm_mat = colMeans(apply(xdiff_array[,,((i_rep-1)*n_reps+1):(i_rep*n_reps)],
                       c(1,2),var))
    sub_mm_mat = as.matrix(1-sub_mm_mat[2:dim(xdiff_array)[2]]/sub_mm_mat[1])
    sub_mm_df_tmp = cbind(known_facs_comb[i,],t(sub_mm_mat))
    sub_mm_df_tmp = melt(sub_mm_df_tmp,
                         id.vars = colnames(sub_mm_df_tmp)[1:3],
                         variable.name = 'scheme')
    mm_df_tmp = rbind(mm_df_tmp,sub_mm_df_tmp)
    
    # time
    sub_mm_mat_t = apply(time_array[,((i_rep-1)*n_reps+1):(i_rep*n_reps)],
                                1,mean)
    sub_mm_mat_t = as.matrix(sub_mm_mat_t[2:dim(xdiff_array)[2]])
    sub_mm_df_tmp_t = cbind(known_facs_comb[i,],t(sub_mm_mat_t))
    sub_mm_df_tmp_t = melt(sub_mm_df_tmp_t,
                         id.vars = colnames(sub_mm_df_tmp_t)[1:3],
                         variable.name = 'scheme')
    mm_df_tmp_t = rbind(mm_df_tmp_t,sub_mm_df_tmp_t)
    
    # iter
    sub_mm_mat_it = apply(iter_array[,((i_rep-1)*n_reps+1):(i_rep*n_reps)],
                          1,mean)
    sub_mm_mat_it = as.matrix(sub_mm_mat_it)
    sub_mm_df_tmp_it = cbind(known_facs_comb[i,],t(sub_mm_mat_it))
    sub_mm_df_tmp_it = melt(sub_mm_df_tmp_it,
                           id.vars = colnames(sub_mm_df_tmp_it)[1:3],
                           variable.name = 'scheme')
    mm_df_tmp_it = rbind(mm_df_tmp_it,sub_mm_df_tmp_it)
    
  }
  xdiff_df = rbind(xdiff_df,mm_df_tmp)
  time_df = rbind(time_df,mm_df_tmp_t)
  iter_df = rbind(iter_df,mm_df_tmp_it)
}

# ANOVA -------------------------------------------------------------------

# covariate
cov_df = xdiff_df
for(i in 1:4){
  cov_df[,i] = as.factor(cov_df[,i])
}
cov_df

aov_cov = aov(value~n*d*rho*scheme,
              data=cov_df)

summary_aov_cov = as.data.frame(summary(aov_cov)[[1]])
summary_aov_cov = cbind(Factors=rownames(summary_aov_cov),
                        summary_aov_cov)
rownames(summary_aov_cov) = NULL
summary_aov_cov = summary_aov_cov[order(-summary_aov_cov$`F value`),]
rownames(summary_aov_cov) = NULL
summary_aov_cov

aov_cov_percent = sum(summary_aov_cov[1:3,]$`Sum Sq`)/sum(summary_aov_cov$`Sum Sq`)
aov_cov_percent

save(aov_cov,
     file = paste0(folder,'aov_cov_group_',
                   as.integer(reps/n_reps),'.rdata'))
write.csv(summary_aov_cov,
          file = paste0(folder,'anova_cov_group',
                        as.integer(reps/n_reps),'.csv'),
          row.names = F)

# time
for(i in 1:4){
  time_df[,i] = as.factor(time_df[,i])
}
time_df

aov_time = aov(value~n*d*rho*scheme,
              data=time_df)

summary_aov_time = as.data.frame(summary(aov_time)[[1]])
summary_aov_time = cbind(Factors=rownames(summary_aov_time),
                        summary_aov_time)
rownames(summary_aov_time) = NULL
summary_aov_time = summary_aov_time[order(-summary_aov_time$`F value`),]
rownames(summary_aov_time) = NULL
summary_aov_time

aov_time_percent = sum(summary_aov_time[1:6,]$`Sum Sq`)/sum(summary_aov_time$`Sum Sq`)
aov_time_percent

save(aov_time,
     file = paste0(folder,'aov_time_group_',
                   as.integer(reps/n_reps),'.rdata'))
write.csv(summary_aov_time,
          file = paste0(folder,'anova_time_group',
                        as.integer(reps/n_reps),'.csv'),
          row.names = F)


# iter
for(i in 1:4){
  iter_df[,i] = as.factor(iter_df[,i])
}
iter_df

aov_iter = aov(value~n*d*rho*scheme,
               data=iter_df)

summary_aov_iter = as.data.frame(summary(aov_iter)[[1]])
summary_aov_iter = cbind(Factors=rownames(summary_aov_iter),
                         summary_aov_iter)
rownames(summary_aov_iter) = NULL
summary_aov_iter = summary_aov_iter[order(-summary_aov_iter$`F value`),]
rownames(summary_aov_iter) = NULL
summary_aov_iter

aov_iter_percent = sum(summary_aov_iter[1:7,]$`Sum Sq`)/sum(summary_aov_iter$`Sum Sq`)
aov_iter_percent

save(aov_iter,
     file = paste0(folder,'aov_iter_group_',
                   as.integer(reps/n_reps),'.rdata'))
write.csv(summary_aov_iter,
          file = paste0(folder,'anova_iter_group',
                        as.integer(reps/n_reps),'.csv'),
          row.names = F)


# summary table -----------------------------------------------------------

# covariate
# three-way tables
agg_cov_df = aggregate(value~d+rho+scheme,data = cov_df, mean)
agg_cov_df
rfmt_cov_df = dcast(data=agg_cov_df,scheme~d+rho)
rfmt_cov_df

#agg_pcadim_df = aggregate(pca_dim~d+rho,data = pca_dim_df, mean)
#rfmt_pcadim_df = dcast(data=agg_pcadim_df,.~d+rho)
#rfmt_pcadim_df = cbind(rfmt_pcadim_df[1,])
#colnames(rfmt_pcadim_df)[1] = colnames(rfmt_cov_df)[1]

rfmt_cov_all_df = rbind(rfmt_cov_df)
rfmt_cov_all_df

write.csv(rfmt_cov_all_df,
          file = paste0(folder,'table_cov.csv'),
          row.names = F)

# time
agg_time_df = aggregate(value~n+d+scheme,data = time_df, mean)
agg_time_df
rfmt_time_df = dcast(data=agg_time_df,scheme~n+d)
rfmt_time_df

#agg_pcadim_df = aggregate(pca_dim~n+d,data = pca_dim_df, mean)
#rfmt_pcadim_df = dcast(data=agg_pcadim_df,.~n+d)
#rfmt_pcadim_df = cbind(rfmt_pcadim_df[1,])
#colnames(rfmt_pcadim_df)[1] = colnames(rfmt_time_df)[1]

rfmt_time_all_df = rbind(rfmt_time_df)
rfmt_time_all_df
write.csv(rfmt_time_all_df,
          file = paste0(folder,'table_time.csv'),
          row.names = F)

# iteration number
iter_df_summary = subset(iter_df_summary,
                 (n%in%n_seq)&(d%in%d_seq)&(rho%in%rho_seq))
iter_melt_df = melt(iter_df_summary,
     id.vars = colnames(iter_df_summary)[1:3],
     variable.name = 'scheme',
     value.name = 'value')
agg_iter_df = aggregate(value~n+d+scheme,data = iter_melt_df, mean)
rfmt_iter_df = dcast(data=agg_iter_df,scheme~n+d)
write.csv(rfmt_iter_df,
          file = paste0(folder,'table_iter.csv'),
          row.names = F)

# combine iter & time data.frame
rfmt_iter_df1 = rfmt_iter_df
for(i in 2:ncol(rfmt_iter_df)){
  rfmt_iter_df1[,i] = as.integer(rfmt_iter_df1[,i])
}

rfmt_iter_time_df = data.frame()
rfmt_iter_time_df = rbind(rfmt_iter_time_df,
                          rfmt_time_summary_df[1,1:17])
for(i in 1:12){
  rfmt_iter_time_df = rbind(rfmt_iter_time_df,
                            rfmt_time_summary_df[i+1,1:17],
                            rfmt_iter_df1[i,1:17])
}

write.csv(rfmt_iter_time_df,
          file = paste0(folder,'table_iter_time.csv'),
          row.names = F)
