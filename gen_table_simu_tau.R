rm(list=ls())
library(reshape2)
library(forecast)

# this script generates the ANOVA table and other summary tables

folder = './save/06-01Simu/'

load(paste0(folder,'summary_tb.rdata'))

# tau ---------------------------------------------------------------------
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

# Tau
# Load Data  ------------------------------------------------------------------
known_facs_comb = expand.grid(n=n_seq,d=d_seq,rho=rho_seq)
unknown_facs_comb = expand.grid(std_err=std_err_seq,
                                betatype=betatype_seq,
                                ytype=ytype_seq)
diff_df = c()

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
  
  # tau
  tau_array = para_comput$tau_mat
  diff_array = (tau_array-tau)^2
  
  mm_df_tmp = c()
  for(i_rep in 1:floor(reps/n_reps)){
    sub_mm_mat = apply(diff_array[,,((i_rep-1)*n_reps+1):(i_rep*n_reps)],
                       c(1,2),mean)
    sub_mm_mat = 1-sub_mm_mat[,2:dim(diff_array)[2]]/sub_mm_mat[,1]
    sub_mm_df_tmp = cbind(known_facs_comb[i,],
                          unknown_facs_comb,
                          sub_mm_mat)
    sub_mm_df_tmp = melt(sub_mm_df_tmp,
                         id.vars = colnames(sub_mm_df_tmp)[1:6],
                         variable.name = 'scheme')
    mm_df_tmp = rbind(mm_df_tmp,sub_mm_df_tmp)
  }
  diff_df = rbind(diff_df,mm_df_tmp)
}

# ANOVA -------------------------------------------------------------------
tau_df = diff_df
for(i in 1:7){
  tau_df[,i] = as.factor(tau_df[,i])
}
tau_df

aov_tau = aov(value~n*d*rho*scheme*betatype*ytype*std_err,
              data=tau_df)

# check anova
# library(car)
# qqPlot(aov_tau$residuals)
# tmp = aggregate(value~n+d+rho+scheme+betatype+ytype+std_err,tau_df,sd)['value']
# c(min(tmp),max(tmp))
# hist(aov_tau$residuals)

summary_aov_tau = as.data.frame(summary(aov_tau)[[1]])
summary_aov_tau = cbind(Factors=rownames(summary_aov_tau),
                        summary_aov_tau)
rownames(summary_aov_tau) = NULL
summary_aov_tau_full = summary_aov_tau[order(-summary_aov_tau$`F value`),]
if(nrow(summary_aov_tau)>20){
  summary_aov_tau = summary_aov_tau_full[1:20,]
}
rownames(summary_aov_tau) = NULL
summary_aov_tau = rbind(summary_aov_tau,summary_aov_tau_full[nrow(summary_aov_tau_full),])
summary_aov_tau

# percent of variance for the top 4 factors
aov_tau_percent = sum(summary_aov_tau_full[1:4,]$`Sum Sq`)/sum(summary_aov_tau_full$`Sum Sq`)
aov_tau_percent

save(aov_tau,
     file = paste0(folder,'aov_mse_group_',
                   as.integer(reps/n_reps),'.rdata'))
write.csv(summary_aov_tau,
          file = paste0(folder,'anova_tau_mse_group',
                        as.integer(reps/n_reps),'.csv'),
          row.names = F)

# summary table -----------------------------------------------------------
# two-way tables: d 
agg_tau_d_df = aggregate(value~d+scheme,data=tau_df,mean)
rfmt_tau_d_df = dcast(data=agg_tau_d_df,scheme~d)
rfmt_tau_d_df = rbind(rfmt_tau_d_df,
                      c('NA',
                        colMeans(rfmt_tau_d_df[,2:5])))
rfmt_tau_d_df

write.csv(rfmt_tau_d_df,
          file = paste0(folder,'table_tau_d.csv'),
          row.names = F)

# two-way tables: ytype
agg_tau_ytype_df = aggregate(value~ytype+scheme,data=tau_df,mean)
rfmt_tau_ytype_df = dcast(data=agg_tau_ytype_df,scheme~ytype)
rfmt_tau_ytype_df = rbind(rfmt_tau_ytype_df,
                          c('NA',colMeans(rfmt_tau_ytype_df[,2:3])))
rfmt_tau_ytype_df
write.csv(rfmt_tau_ytype_df,
          file = paste0(folder,'table_tau_ytype.csv'),
          row.names = F)

# two-way tables: rho
agg_tau_rho_df = aggregate(value~rho+scheme,data=tau_df,mean)
rfmt_tau_rho_df = dcast(data=agg_tau_rho_df,scheme~rho)
rfmt_tau_rho_df = rbind(rfmt_tau_rho_df,
                          c('NA',colMeans(rfmt_tau_rho_df[,2:4])))
rfmt_tau_rho_df
write.csv(rfmt_tau_rho_df,
          file = paste0(folder,'table_tau_rho.csv'),
          row.names = F)

# two-way tables: n
agg_tau_n_df = aggregate(value~n+scheme,data=tau_df,mean)
rfmt_tau_n_df = dcast(data=agg_tau_n_df,scheme~n)
rfmt_tau_n_df = rbind(rfmt_tau_n_df,
                      c('NA',colMeans(rfmt_tau_n_df[,2:5])))
rfmt_tau_n_df
write.csv(rfmt_tau_n_df,
          file = paste0(folder,'table_tau_n.csv'),
          row.names = F)

# four-way tables
agg_tau_df = aggregate(value~d+rho+ytype+scheme,data = tau_df, mean)
agg_tau_df
rfmt_tau_df = dcast(data=agg_tau_df,ytype+scheme~d+rho)
rfmt_tau_df

rfmt_tau_all_df = rbind(rfmt_tau_df)
rfmt_tau_all_df
write.csv(rfmt_tau_all_df,
          file = paste0(folder,'table_tau.csv'),
          row.names = F)

# generate the aggregate table for pca dimension k
pca_dim_df = subset(pca_dim_df,
       (n%in%n_seq)&(d%in%d_seq)&(rho%in%rho_seq))

pca_dim_df_melt = melt(pca_dim_df,
                       id.vars = colnames(pca_dim_df)[1:3],
                       variable.name = 'scheme')
for(i in 1:4){
  pca_dim_df_melt[,i] = as.factor(pca_dim_df_melt[,i])
}
agg_pcadim_df = aggregate(value~d+rho+scheme,
                          data = pca_dim_df_melt, mean)

rfmt_pcadim_df = dcast(data=agg_pcadim_df,scheme~d+rho)
colnames(rfmt_pcadim_df)[1] = colnames(rfmt_tau_df)[1]
rfmt_pcadim_df

write.csv(rfmt_pcadim_df,
          file = paste0(folder,'table_pcadim.csv'),
          row.names = F)

library(ggplot2)
library(scales)

dim_names <- c(
  '10' = "d=10",
  '50' = "d=50",
  '90' = "d=90",
  '180' = "d=180"
)

library('latex2exp')
barplt = ggplot(agg_pcadim_df, 
                aes(x=rho, 
                    y=value,
                    fill=scheme)) +
  geom_bar(position = "dodge",
           stat = 'identity') +
  facet_wrap(~d, 
             scales = 'free_y',
             labeller = as_labeller(dim_names),
             nrow=1)+
  scale_fill_discrete(name='',
                      labels=unname(TeX(c('$\\gamma_k=0.5$',
                                          '$\\gamma_k=0.7$',
                                          '$\\gamma_k=0.9$',
                                          'Kaiser'))))+
  xlab(TeX('$\\rho$'))+
  ylab(TeX('$k$')) + 
  theme_bw()

ggsave(file.path(paste0(folder,'pcadim_barplot.pdf')),
       plot=barplt,width = 8,height = 4)
barplt

