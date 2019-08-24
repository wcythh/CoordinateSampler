rm(list=ls(all=TRUE))
install.packages("mcmcse")
library(MASS)
library(mcmcse)
library(parallel)
d <- 2
kappa <- 20
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

ZZ_Coef <- function(x,v,k)
{
  if(k == 1)
  {
    return(c(4*kappa*v[1]^4, 12*kappa*x[1]*v[1]^3 -4*kappa*v[1]^2*v[2],
             2*v[1]^2+12*kappa*x[1]^2*v[1]^2-4*kappa*x[1]*v[1]*v[2]-4*kappa*x[2]*v[1]^2,
             2*x[1]*v[1]-2*v[1]+4*kappa*v[1]*x[1]^3 - 4*kappa*v[1]*x[1]*x[2]))
  }
  if(k == 2)
  {
    return(c(-2*kappa*v[1]^2*v[2], 2*kappa*v[2]*(v[2]-2*x[1]*v[1]),
             2*kappa*v[2]*(x[2]-x[1]^2)))
  }
}

FirstEventTimeZZ <- function(x,v,k)
{
  Ttemp <- 0
  if(k == 1)
  {
    Stop <- FALSE
    while(!Stop)
    {
      a <- ZZ_Coef((x+Ttemp*v),v,k)
      aPositive <- PositivePart(a)
      aSum <- sum(aPositive)
      u <- -log(runif(1))
      if((u/aSum) <= 1)
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
        t0 <- (4*u/aSum -3)^0.25
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
    Stop <- FALSE
    while(!Stop)
    {
      a <- ZZ_Coef((x+Ttemp*v),v,k)
      aPositive <- PositivePart(a)
      aSum <- sum(aPositive)
      if(aSum == 0)
      {
        Stop <- TRUE
        return(1e7)
      }
      u <- -log(runif(1))
      if((u/aSum) <= 1)
      {
        t0 <- u/aSum
        valTrue <- a[1]*t0^2 + a[2]*t0^1 + a[3]
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
        t0 <- (3*u/aSum -2)^(1/3)
        valTrue <- a[1]*t0^2 + a[2]*t0 + a[3]
        valUpper <- aSum*t0^2
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
}

ZZPDMP <- function(mmm)
{
  set.seed(mmm*1e5)
  Niter_zz <- 5e5
  Xtrace_zz <- matrix(0,nrow=Niter_zz,ncol=d)
  Vtrace_zz <- matrix(0,nrow=Niter_zz,ncol=d)
  Ttrace_zz <- rep(0,Niter_zz)
  Xtrace_zz[1,] <- rnorm(d)
  Vtrace_zz[1,] <- sample(c(-1,1),size=d,replace=TRUE)
  Total_zz <- 0
  time_zz <- proc.time()
  for(i in 2:Niter_zz)
  {
    Ttemp <- rep(0,d)
    for(k in 1:d)
    {
      Ttemp[k] <- FirstEventTimeZZ(Xtrace_zz[(i-1),],Vtrace_zz[(i-1),],k)
    }
    Ind_zz <- which.min(Ttemp)
    Ttrace_zz[i] <- Ttemp[Ind_zz]
    Xtrace_zz[i,] <- Xtrace_zz[(i-1),] + Ttrace_zz[i] * Vtrace_zz[(i-1),]
    Vtrace_zz[i,] <- Vtrace_zz[(i-1),]
    Vtrace_zz[i,Ind_zz] <- -Vtrace_zz[i,Ind_zz]
    #if(i %% (Niter_zz/10) == 0)
    #  print(paste("Iteration = ", toString(i)))
    if(Ind_zz == 1)
      Total_zz <- Total_zz + 1
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
  
  lp_zz <- apply(Xsamp_zz,1,function(x) U(x))

  return(c(ess(Xsamp_zz), ess(lp_zz), time_zz[[1]], ess(Xsamp_zz)/time_zz[[1]], ess(lp_zz)/time_zz[[1]]))
}

num_cores <- detectCores()

cl <- makeCluster(num_cores,type = "FORK")

X_zz <- parSapply(cl, 1:40, function(x) ZZPDMP(x))

save(X_zz, file="X_zz.RData")
stopCluster(cl)
X_zz <- t(X_zz)
print(apply(X_zz,2,mean))
