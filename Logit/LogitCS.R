install.packages("mcmcse")
library(mcmcse)
library(MASS)
N <- 100
d <- 10

R <- get(load(file="R.RData"))
S <- get(load(file="S.RData"))

lambda_0 <- 0.1
C <- apply(abs(R),2,sum) + lambda_0

train <- data.frame(R,S)
model <- glm(S ~.-1,family=binomial(link='logit'),data=train)
model$coefficients

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

ESS_cs <- matrix(0,nrow=50,ncol=d)
ESS_LP_cs <- matrix(0,nrow=50,ncol=d)
Time_cs <- rep(0,50)
for(Iter in 1:50)
{
  Niter_cs <- 5e3*d
  Xtrace_cs <- matrix(0,nrow=Niter_cs,ncol=d)
  Vtrace_cs <- matrix(0,nrow=Niter_cs,ncol=d)
  Xtrace_cs[1,] <- rep(0,d)
  Ind_cs <- sample(1:d,size=1)
  Vtrace_cs[1,Ind_cs] <- sample(c(-1,1),size=1)
  Ttrace_cs <- rep(0,Niter_cs)
  time_cs <- proc.time()
  for(i in 2:Niter_cs)
  {
    Ttrace_cs[i] <- FirstEventTime(Xtrace_cs[(i-1),], Vtrace_cs[(i-1),],Ind_cs)
    Xtrace_cs[i,] <- Xtrace_cs[(i-1),] + Ttrace_cs[i] * Vtrace_cs[(i-1),]
    p <- cbind(0,grad_U(Xtrace_cs[i,]))
    p <- c(apply(p,1,max), apply(-p,1,max)) + lambda_0
    Ind_cs <- sample(c(1:(2*d)),size=1,prob=p)
    if(Ind_cs <= d)
      Vtrace_cs[i,Ind_cs] = -1
    else
    {
      Ind_cs <- Ind_cs - d
      Vtrace_cs[i,(Ind_cs)] = 1
    }
    #if(i %% (Niter_cs/10) ==0)
    #  print(i)
  }
  time_cs <- proc.time() - time_cs
  
  Num_cs <- 1e4
  fact_cs <- sum(Ttrace_cs)/Num_cs
  Tcumsum_cs <- cumsum(Ttrace_cs)
  Xsamp_cs <- matrix(0,nrow=(Num_cs),ncol=d)
  Vsamp_cs <- matrix(0,nrow=(Num_cs),ncol=d)
  for(i in 1:Num_cs)
  {
    m <- max(which(Tcumsum_cs < fact_cs*i))
    delta <- fact_cs*i - Tcumsum_cs[m]
    Xsamp_cs[i,] <- Xtrace_cs[m,] + Vtrace_cs[m,]*delta
    Vsamp_cs[i,] <- Vtrace_cs[m,]
  }
  
  lp_cs <- apply(Xsamp_cs,1,U)
  
  cat(paste("Iter = ", toString(Iter),sep=""),
      file="outputCS.txt",sep="\n",append=TRUE)
  cat(paste("........", "Time-consuming = ", toString(time_cs[[1]]), sep=""),
      file="outputCS.txt",sep="\n",append=TRUE)
  cat(paste("........",
            "MinEss=", toString(min(ess(Xsamp_cs))),
            " MeanEss=",toString(mean(ess(Xsamp_cs))),
            " MedianEss=", toString(median(ess(Xsamp_cs))),
            " MaxEss=", toString(max(ess(Xsamp_cs))),
            sep=""),
      file="outputCS.txt",sep="\n",append=TRUE)
  cat(paste("........",
            "MinEss/s=", toString(min(ess(Xsamp_cs))/time_cs[[1]]),
            " MeanEss/s=",toString(mean(ess(Xsamp_cs))/time_cs[[1]]),
            " MedianEss/s=", toString(median(ess(Xsamp_cs))/time_cs[[1]]),
            " MaxEss/s=", toString(max(ess(Xsamp_cs))/time_cs[[1]]),
            sep=""),
      file="outputCS.txt",sep="\n",append=TRUE)
  cat(paste("........",
            "Ess(lp)=", toString(ess(lp_cs)),
            " Ess(lp)/s=",toString(ess(lp_cs)/time_cs[[1]]),
            sep=""),
      file="outputCS.txt",sep="\n",append=TRUE)
  
  print(c(Iter, time_cs[[1]], ess(lp_cs), ess(lp_cs)/time_cs[[1]]))
  print(c(Iter, apply(Xsamp_cs,2,mean)))
  save(Xsamp_cs,file=paste("CS_",toString(Iter),".RData",sep=""))
  plot(density(Xsamp_cs[,1]))
  curve(dnorm(x),col="red",add=TRUE)
  acf(lp_cs)
  ESS_cs[Iter,] = ess(Xsamp_cs)
  ESS_LP_cs[Iter] = ess(lp_cs)
  Time_cs[Iter] = time_cs[[1]]
}
save(ESS_cs, file="ESS_CS.RData")
save(ESS_LP_cs, file="ESS_LP_CS.RData")
save(Time_cs, file="Time_CS.RData")


