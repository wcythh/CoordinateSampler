Compare_KS_min <- matrix(0, nrow=40, ncol=4)
Compare_KS_mean <- matrix(0, nrow=40, ncol=4)
Compare_KS_median <- matrix(0, nrow=40, ncol=4)
Compare_KS_max <- matrix(0, nrow=40, ncol=4)
Compare_time <- matrix(0, nrow=40, ncol=4)
Compare_grad <- matrix(0,nrow=40,ncol=4)

ZZResult <- get(load(file="Results/ResultZZ.RData"))
CSResult <- get(load(file="Results/ResultCS.RData"))
BPSResult <- get(load(file="Results/ResultBPS.RData"))
HMCResult <- get(load(file="Results/ResultHMC.RData"))

for(j in 1:40)
{
  Compare_KS_min[j,] <- c(min(CSResult[[j]]$Xsamp_ks),min(ZZResult[[j]]$Xsamp_ks), min(BPSResult[[j]]$Xsamp_ks), min(HMCResult[[j]]$Xsamp_ks))
  Compare_KS_mean[j,] <- c(mean(CSResult[[j]]$Xsamp_ks),mean(ZZResult[[j]]$Xsamp_ks), mean(BPSResult[[j]]$Xsamp_ks), mean(HMCResult[[j]]$Xsamp_ks))
  Compare_KS_median[j,] <- c(median(CSResult[[j]]$Xsamp_ks),median(ZZResult[[j]]$Xsamp_ks), median(BPSResult[[j]]$Xsamp_ks), median(HMCResult[[j]]$Xsamp_ks))
  Compare_KS_max[j,] <- c(max(CSResult[[j]]$Xsamp_ks),max(ZZResult[[j]]$Xsamp_ks), max(BPSResult[[j]]$Xsamp_ks), max(HMCResult[[j]]$Xsamp_ks))
  Compare_time[j,] <- c(CSResult[[j]]$time[[1]], ZZResult[[j]]$time[[1]], BPSResult[[j]]$time[[1]], HMCResult[[j]]$time[[1]])
  Compare_grad[j,] <- c(CSResult[[j]]$grad_num, ZZResult[[j]]$grad_num, BPSResult[[j]]$grad_num, HMCResult[[j]]$grad_num)
}

boxplot(Compare_KS_min)
boxplot(Compare_KS_mean)
boxplot(Compare_KS_median)
boxplot(Compare_KS_max)
boxplot(Compare_time)
boxplot(Compare_grad)

apply(Compare_KS_min,2,mean)
apply(Compare_KS_mean,2,mean)
apply(Compare_KS_median,2,mean)
apply(Compare_KS_max,2,mean)


apply(Compare_time,2,mean)
apply(Compare_grad,2,mean)








