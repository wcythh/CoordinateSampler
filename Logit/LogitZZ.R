install.packages("mcmcse")
library(mcmcse)
library(MASS)
N <- 100
d <- 10

R <- get(load(file="R.RData"))
S <- get(load(file="S.RData"))

lambda_0 <- 0.1
C <- apply(abs(R),2,sum) + lambda_0

#train <- data.frame(R,S)
#model <- glm(S ~.-1,family=binomial(link='logit'),data=train)
#model$coefficients

grad_U <- function(x)
{
  val <- as.vector(R %*% x)
  val <- 1.0/(1+exp(-val)) - S
  val <- as.vector(t(R) %*% val)
  return(val)
}

U <- function(x)
{
  val <- as.vector(R%*%x)
  val <- S*val - log(1+exp(val))
  return(sum(val))
}

FirstEventTime <- function(x,v,k)
{
  Ttemp = 0.0
  Stop = FALSE
  while(!Stop)
  {
    t0 = -log(runif(1))/C[k]
    Ttemp = Ttemp + t0
    val = v[k] * (grad_U(x+Ttemp*v))[k]
    val <- max(val,0) + lambda_0
    if(val >= 0)
    {
      u <- runif(1)
      if(u < ((val)/C[k]))
      {
        Stop = TRUE
      }
    }
  }
  return(Ttemp)
}
ESS_zz <- matrix(0,nrow=50,ncol=d)
ESS_LP_zz <- matrix(0,nrow=50,ncol=d)
Time_zz <- rep(0,50)

for(Iter in 1:50)
{
  Niter_zz <- 1e4
  Xtrace_zz <- matrix(0,nrow=Niter_zz,ncol=d)
  Vtrace_zz <- matrix(0,nrow=Niter_zz,ncol=d)
  Ttrace_zz <- rep(0,Niter_zz)
  Xtrace_zz[1,] <- rep(0,d)
  #Xtrace_zz[1,] <- model$coefficients
  Vtrace_zz[1,] <- sample(c(-1,1),size=d,replace=TRUE)
  
  time_zz <- proc.time()
  for(i in 2:Niter_zz)
  {
    Tt <- rep(0,d)
    for(k in 1:d)
    {
      Tt[k] <- FirstEventTime(Xtrace_zz[(i-1),], Vtrace_zz[(i-1),], k)
    }
    Ind = which.min(Tt)
    Xtrace_zz[i,] <- Xtrace_zz[(i-1),] + Tt[Ind] * Vtrace_zz[(i-1),]
    Vtrace_zz[i,] <- Vtrace_zz[(i-1),]
    Vtrace_zz[i,Ind] <- -Vtrace_zz[i,Ind]
    Ttrace_zz[i] <- Tt[Ind]
    #if(i %% (Niter_zz/10) == 0)
    #  print(i)
  }
  time_zz <- proc.time() - time_zz
  
  Num_zz <- 1e4
  fact_zz <- sum(Ttrace_zz)/Num_zz
  Tcumsum_zz <- cumsum(Ttrace_zz)
  Xsamp_zz <- matrix(0,nrow=Num_zz,ncol=d)
  Vsamp_zz <- matrix(0,nrow=Num_zz,ncol=d)
  for(i in 1:Num_zz)
  {
    m <- max(which(Tcumsum_zz < i*fact_zz))
    delta <- i*fact_zz - Tcumsum_zz[m]
    Xsamp_zz[i,] <- Xtrace_zz[m,] + delta*Vtrace_zz[m,]
    Vsamp_zz[i,] <- Vtrace_zz[m,]
  }
  lp_zz <- apply(Xsamp_zz,1,U)
  
  cat(paste("Iter = ", toString(Iter),sep=""),
      file="outputZZ.txt",sep="\n",append=TRUE)
  cat(paste("........", "Time-consuming = ", toString(time_zz[[1]]), sep=""),
      file="outputZZ.txt",sep="\n",append=TRUE)
  cat(paste("........",
            "MinEss=", toString(min(ess(Xsamp_zz))),
            " MeanEss=",toString(mean(ess(Xsamp_zz))),
            " MedianEss=", toString(median(ess(Xsamp_zz))),
            " MaxEss=", toString(max(ess(Xsamp_zz))),
            sep=""),
      file="outputZZ.txt",sep="\n",append=TRUE)
  cat(paste("........",
            "MinEss/s=", toString(min(ess(Xsamp_zz))/time_zz[[1]]),
            " MeanEss/s=",toString(mean(ess(Xsamp_zz))/time_zz[[1]]),
            " MedianEss/s=", toString(median(ess(Xsamp_zz))/time_zz[[1]]),
            " MaxEss/s=", toString(max(ess(Xsamp_zz))/time_zz[[1]]),
            sep=""),
      file="outputZZ.txt",sep="\n",append=TRUE)
  cat(paste("........",
            "Ess(lp)=", toString(ess(lp_zz)),
            " Ess(lp)/s=",toString(ess(lp_zz)/time_zz[[1]]),
            sep=""),
      file="outputZZ.txt",sep="\n",append=TRUE)
  
  print(c(Iter, time_zz[[1]], ess(lp_zz), ess(lp_zz)/time_zz[[1]]))
  print(c(Iter, apply(Xsamp_zz,2,mean)))
  save(Xsamp_zz,file=paste("ZZ_",toString(Iter),".RData",sep=""))
  plot(density(Xsamp_zz[,1]))
  curve(dnorm(x),col="red",add=TRUE)
  acf(lp_zz)
  ESS_zz[Iter,] = ess(Xsamp_zz)
  ESS_LP_zz[Iter] = ess(lp_zz)
  Time_zz[Iter] = time_zz[[1]]
}
save(ESS_zz, file="ESS_ZZ.RData")
save(ESS_LP_zz, file="ESS_LP_ZZ.RData")
save(Time_zz, file="Time_ZZ.RData")



