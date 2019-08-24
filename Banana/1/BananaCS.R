rm(list = ls(all=TRUE))
install.packages("mcmcse")
library(parallel)
library(MASS)
library(mcmcse)
kappa <- 2^{0}
U <- function(x)
{
  return((x[1]-1)^2+kappa*(x[2]-x[1]^2)^2)
}
PositivePart <- function(a)
{
  return(apply(rbind(a,rep(0,length(a))),2,max))
}

grad_U <- function(x)
{
  return(c(2*(x[1]-1)+4*kappa*x[1]*(x[1]^2-x[2]), 2*kappa*(x[2] - x[1]^2)))
}

CS_Coef <- function(x,v,k)
{
  if(k == 1)
  {
    return(c(4*kappa*v[1]^4,12*kappa*x[1]*v[1]^3,
             2*v[1]^2 + 12*kappa*x[1]^2*v[1]^2 - 4*kappa*x[2]*v[1]^2,
             2*x[1]*v[1]-2*v[1]+4*kappa*v[1]*x[1]^3 - 4*kappa*x[1]*x[2]*v[1]))
  }
  if(k == 2)
  {
    return(c(2*kappa*v[2]^2, 2*kappa*v[2]*(x[2]-x[1]^2)))
  }
}

FirstEventTimeCS <- function(x,v,k)
{
  Ttemp <- 0
  if(k == 1)
  {
    Stop <- FALSE
    while(!Stop)
    {
      a <- CS_Coef((x+Ttemp*v), v, k)
      aPositive <- PositivePart(a)
      aSum <- sum(aPositive)
      u <- -log(runif(1))
      if(u/aSum <= 1)
      {
        t0 <- u/aSum
        valTrue <- a[1]*t0^3 + a[2]*t0^2 + a[3]*t0 + a[4]
        valUpper <- aSum
        p <- runif(1)
        if(p < (valTrue/valUpper))
        {
          Ttemp <- Ttemp + t0
          Stop <- TRUE
          return(Ttemp)
        }
        else
        {
          Ttemp <- Ttemp + t0
        }
      }
      else
      {
        t0 <- (4*u/aSum - 3)^0.25
        valTrue <- a[1]*t0^3 + a[2]*t0^2 + a[3]*t0 + a[4]
        valUpper <- aSum*t0^3
        p <- runif(1)
        if(p < (valTrue/valUpper))
        {
          Ttemp <- Ttemp + t0
          Stop <- TRUE
          return(Ttemp)
        }
        else
        {
          Ttemp <- Ttemp + t0
        }
      }
      
    }
  }
  if(k == 2)
  {
    a = CS_Coef(x,v,k)
    u = -log(runif(1))
    if(a[2] > 0)
      return(sqrt((2*u*a[1]+a[2]^2))/a[1] - a[2]/a[1])
    else
      return(sqrt(2*u/a[1]) - a[2]/a[1])
  }
}

d <- 2
ESS_CS <- matrix(0,nrow=10,ncol=5)
CSPDMP <- function(mmm)
{
  set.seed(mmm*1e5)
  Niter_cs <- 1e5
  Xtrace_cs <- matrix(0,nrow=Niter_cs,ncol=d)
  Vtrace_cs <- matrix(0,nrow=Niter_cs,ncol=d)
  Ttrace_cs <- rep(0,Niter_cs)
  
  Xtrace_cs[1,] <- rnorm(d)
  Ind_cs <- sample(1:d,size=1)
  Vtrace_cs[1,Ind_cs] <- sample(c(-1,1),size=1)
  Total_cs <- 0
  time_cs <- proc.time()
  for(i in 2:Niter_cs)
  {
    Ttrace_cs[i] <- FirstEventTimeCS(Xtrace_cs[(i-1),], Vtrace_cs[(i-1),], Ind_cs)
    Xtrace_cs[i,] <- Xtrace_cs[(i-1),] + Ttrace_cs[i] * Vtrace_cs[(i-1),]
    p <- cbind(0,grad_U(Xtrace_cs[i,]))
    p <- c(apply(p,1,max), apply(-p,1,max))
    Ind_cs <- sample(c(1:(2*d)),size=1,prob=p)
    if(Ind_cs <= d)
      Vtrace_cs[i,Ind_cs] = -1
    else
    {
      Ind_cs <- Ind_cs - d
      Vtrace_cs[i,(Ind_cs)] = 1
    }
    #if(i %% (Niter_cs/10) == 0)
    #  print(paste("Iteration = ", toString(i)))
    if(Ind_cs == 1)
      Total_cs <- Total_cs + 1
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
  
  lp_cs <- apply(Xsamp_cs,1,function(x) U(x))

  return(c(ess(Xsamp_cs), ess(lp_cs), time_cs[[1]], ess(Xsamp_cs)/time_cs[[1]], ess(lp_cs)/time_cs[[1]]))
}

num_cores <- detectCores()

cl <- makeCluster(num_cores,type = "FORK")

X_cs <- parSapply(cl, 1:40, function(x) CSPDMP(x))

save(X_cs, file="X_cs.RData")
stopCluster(cl)
X_cs <- t(X_cs)
print(apply(X_cs,2,mean))
