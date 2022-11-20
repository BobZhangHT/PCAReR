rm(list=ls())
library(reshape2)
library(forecast)

# anova table for the latest factorial design
folder = './save/11-29V1/'

# Covariate ---------------------------------------------------------------
load(paste0(folder,'summary_tb.rdata'))
# analyze the known factors (rho & zeta & n & scheme) for covariates

cov_df = melt(diff_ratio_df,id.vars = c('n','rnd','rho'),
     variable.name = 'scheme')
cov_df[,1] = as.factor(cov_df[,1])
cov_df[,2] = as.factor(cov_df[,2])
cov_df[,3] = as.factor(cov_df[,3])
cov_df[,4] = as.factor(cov_df[,4])
cov_df[,5] = 1-cov_df[,5]

aov_cov = aov(value~n*rnd*rho*scheme,
              data=cov_df)
summary_aov_cov = as.data.frame(summary(aov_cov)[[1]])
summary_aov_cov = cbind(Factors=rownames(summary_aov_cov),
                            summary_aov_cov)
rownames(summary_aov_cov) = NULL
summary_aov_cov = summary_aov_cov[order(-summary_aov_cov$`Mean Sq`),]
if(nrow(summary_aov_cov)>20){
  summary_aov_cov = summary_aov_cov[1:20,]
}
rownames(summary_aov_cov) = NULL
summary_aov_cov
write.csv(summary_aov_cov,
          file = paste0(folder,'anova_cov.csv'),
          row.names = F)

# Time --------------------------------------------------------------------
time_aov_df = melt(time_df,id.vars = c('n','rnd','rho'),
               variable.name = 'scheme')
time_aov_df[,1] = as.factor(time_aov_df[,1])
time_aov_df[,2] = as.factor(time_aov_df[,2])
time_aov_df[,3] = as.factor(time_aov_df[,3])
time_aov_df[,4] = as.factor(time_aov_df[,4])  
aov_time = aov(value~n*rnd*rho*scheme,
               data=time_aov_df)
summary_aov_time = as.data.frame(summary(aov_time)[[1]])
summary_aov_time = cbind(Factors=rownames(summary_aov_time),
                            summary_aov_time)
rownames(summary_aov_time) = NULL
summary_aov_time = summary_aov_time[order(-summary_aov_time$`Mean Sq`),]
rownames(summary_aov_time) = NULL
summary_aov_time
write.csv(summary_aov_time,
          file = paste0(folder,'anova_time.csv'),
          row.names = F)

# Tau --------------------------------------------------------------------
tau_df = melt(mse_df,id.vars = colnames(mse_df)[1:6],
                   variable.name = 'scheme')
for(i in 1:7){
  tau_df[,i] = as.factor(tau_df[,i])
}
tau_df[,8] = 1-tau_df[,8]
aov_tau = aov(value~n*rnd*rho*scheme*betatype*ytype*std_err,
               data=tau_df)
summary_aov_tau = as.data.frame(summary(aov_tau)[[1]])
summary_aov_tau = cbind(Factors=rownames(summary_aov_tau),
                        summary_aov_tau)
rownames(summary_aov_tau) = NULL
summary_aov_tau = summary_aov_tau[order(-summary_aov_tau$`Mean Sq`),]
if(nrow(summary_aov_tau)>20){
  summary_aov_tau = summary_aov_tau[1:20,]
}
rownames(summary_aov_tau) = NULL
summary_aov_tau
write.csv(summary_aov_tau,
          file = paste0(folder,'anova_tau.csv'),
          row.names = F)
