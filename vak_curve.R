# generate the vak plot
library(latex2exp) 

d = 200
k = seq(1,d)
pa = 0.05
ak = qchisq(pa,k)
vak = pchisq(ak,df=k+2)/pchisq(ak,df=k)

pdf('./save/vak.pdf',
    height = 3,width = 4)
par(mar=c(4.5,4.5,2,2))
plot(k,vak,type='l',
     ylab=TeX('v_{a_k}'),
     xlab = 'k')
dev.off()