Len <- c(2^(-2), 2^(-1), 2^0, 2^1, 2^2, 2^3,12, 2^4,20,24,26,2^5)
length(Len)
val1 <- rep(0,12)
val2 <- rep(0,12)
val3 <- rep(0,12)



Result1_cs <- t(get(load("Banana/0_25/X_cs.RData")))
Result1_zz <- t(get(load("Banana/0_25/X_zz.RData")))
val1[1] <- mean(Result1_cs[,5])/mean(Result1_zz[,5])
val2[1] <- mean(Result1_cs[,6])/mean(Result1_zz[,6])
val3[1] <- mean(Result1_cs[,7])/mean(Result1_zz[,7])


Result2_cs <- t(get(load("Banana/0_5/X_cs.RData")))
Result2_zz <- t(get(load("Banana/0_5/X_zz.RData")))
val1[2] <- mean(Result2_cs[,5])/mean(Result2_zz[,5])
val2[2] <- mean(Result2_cs[,6])/mean(Result2_zz[,6])
val3[2] <- mean(Result2_cs[,7])/mean(Result2_zz[,7])


Result3_cs <- t(get(load("Banana/1/X_cs.RData")))
Result3_zz <- t(get(load("Banana/1/X_zz.RData")))
val1[3] <- mean(Result3_cs[,5])/mean(Result3_zz[,5])
val2[3] <- mean(Result3_cs[,6])/mean(Result3_zz[,6])
val3[3] <- mean(Result3_cs[,7])/mean(Result3_zz[,7])


Result4_cs <- t(get(load("Banana/2/X_cs.RData")))
Result4_zz <- t(get(load("Banana/2/X_zz.RData")))
val1[4] <- mean(Result4_cs[,5])/mean(Result4_zz[,5])
val2[4] <- mean(Result4_cs[,6])/mean(Result4_zz[,6])
val3[4] <- mean(Result4_cs[,7])/mean(Result4_zz[,7])


Result5_cs <- t(get(load("Banana/4/X_cs.RData")))
Result5_zz <- t(get(load("Banana/4/X_zz.RData")))
val1[5] <- mean(Result5_cs[,5])/mean(Result5_zz[,5])
val2[5] <- mean(Result5_cs[,6])/mean(Result5_zz[,6])
val3[5] <- mean(Result5_cs[,7])/mean(Result5_zz[,7])

Result6_cs <- t(get(load("Banana/8/X_cs.RData")))
Result6_zz <- t(get(load("Banana/8/X_zz.RData")))
val1[6] <- mean(Result6_cs[,5])/mean(Result6_zz[,5])
val2[6] <- mean(Result6_cs[,6])/mean(Result6_zz[,6])
val3[6] <- mean(Result6_cs[,7])/mean(Result6_zz[,7])

Result7_cs <- t(get(load("Banana/12/X_cs.RData")))
Result7_zz <- t(get(load("Banana/12/X_zz.RData")))
val1[7] <- mean(Result7_cs[,5])/mean(Result7_zz[,5])
val2[7] <- mean(Result7_cs[,6])/mean(Result7_zz[,6])
val3[7] <- mean(Result7_cs[,7])/mean(Result7_zz[,7])

Result8_cs <- t(get(load("Banana/16/X_cs.RData")))
Result8_zz <- t(get(load("Banana/16/X_zz.RData")))
val1[8] <- mean(Result8_cs[,5])/mean(Result8_zz[,5])
val2[8] <- mean(Result8_cs[,6])/mean(Result8_zz[,6])
val3[8] <- mean(Result8_cs[,7])/mean(Result8_zz[,7])

Result9_cs <- t(get(load("Banana/20/X_cs.RData")))
Result9_zz <- t(get(load("Banana/20/X_zz.RData")))
val1[9] <- mean(Result9_cs[,5])/mean(Result9_zz[,5])
val2[9] <- mean(Result9_cs[,6])/mean(Result9_zz[,6])
val3[9] <- mean(Result9_cs[,7])/mean(Result9_zz[,7])

Result10_cs <- t(get(load("Banana/24/X_cs.RData")))
Result10_zz <- t(get(load("Banana/24/X_zz.RData")))
val1[10] <- mean(Result10_cs[,5])/mean(Result10_zz[,5])
val2[10] <- mean(Result10_cs[,6])/mean(Result10_zz[,6])
val3[10] <- mean(Result10_cs[,7])/mean(Result10_zz[,7])

Result11_cs <- t(get(load("Banana/26/X_cs.RData")))
Result11_zz <- t(get(load("Banana/26/X_zz.RData")))
val1[11] <- mean(Result11_cs[,5])/mean(Result11_zz[,5])
val2[11] <- mean(Result11_cs[,6])/mean(Result11_zz[,6])
val3[11] <- mean(Result11_cs[,7])/mean(Result11_zz[,7])


Result12_cs <- t(get(load("Banana/32/X_cs.RData")))
Result12_zz <- t(get(load("Banana/32/X_zz.RData")))
val1[12] <- mean(Result12_cs[,5])/mean(Result12_zz[,5])
val2[12] <- mean(Result12_cs[,6])/mean(Result12_zz[,6])
val3[12] <- mean(Result12_cs[,7])/mean(Result12_zz[,7])



Ind <- c(1:12)
plot(log(Len[Ind],2),val1[Ind],pch=20,type="o",col="red",ylim=c(0,3),
     xlab = expression(log[2](kappa)), ylab="Ratio of ESS per second")
lines(log(Len[Ind],2),val2[Ind],pch=20,type="o",col="blue")
lines(log(Len[Ind],2),val3[Ind],pch=20,type="o",col="green")
legend(-2.27,3.11, col=c("red","blue", "green"),
       legend = c("first component","second component","log-likelihood"),lty=1)
