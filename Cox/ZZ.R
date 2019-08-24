library(MASS)
library(mcmcse)
set.seed(123456)
d0 <- 20
s0 <- 1/(d0^2)
beta0 <- 1/6
sigma2 <- 1.91
mu0 <- log(126) - 0.5*sigma2

# k = d*(i-1) + j
Sigma0 <- matrix(0,nrow=(d0^2),ncol=(d0^2))
for(k in 1:(d0^2))
{
  for(l in 1:(d0^2))
  {
    val1 <- c((k %/% d0) + 1, (k %% d0))
    val2 <- c((l %/% d0) + 1, (l %% d0))
    if(val1[2] == 0)
      val1[2] <- d0
    if(val2[2] == 0)
      val2[2] <- d0
    delta <- sqrt(sum((val1-val2)^2))
    Sigma0[k,l] <- sigma2 * exp(-delta/(d0*beta0))
  }
}

Sigma_Inv <- solve(Sigma0)
X0 <- mvrnorm(1,mu=rep(mu0,(d0^2)),Sigma=Sigma0)

par(mar=c(2,2,2,2))
plot(X0,type="o",pch=20,lwd=2)

Lambda0 <- exp(X0)
Y0 <- rep(0,(d0^2))
for(k in 1:(d0^2))
{
  Y0[k] <- rpois(1,lambda=(s0*Lambda0[k]))
}
plot(Y0,type="o",pch=20,lwd=2)

d <- d0^2
lambda0 <- 0.1

U <- function(theta)
{
  val <- -sum(Y0 * theta) + s0*sum(exp(theta)) + 
    0.5*as.numeric(t(theta-mu0)%*%Sigma_Inv%*%(theta-mu0))
  return(val)
}

grad_U <- function(theta)
{
  val <- -Y0 + s0*exp(theta) + as.vector(Sigma_Inv%*%(theta-mu0))
  return(val)
}

FirstEventTime1 <- function(a,b)
{
  u <- -log(runif(1))
  if(a > 0)
  {
    if(b >0)
    {
      val <- sqrt(b^2+2*a*u)/a - b/a
      return(val)
    }
    else
    {
      val <- sqrt(2*u/a) - b/a
      return(val)
    }
  }
  if(a == 0)
  {
    if(b > 0)
    {
      val <- u/b
      return(val)
    }
    else
    {
      return(Inf)
    }
  }
  if(a < 0)
  {
    if(b <= 0)
    {
      val <- Inf
      return(val)
    }
    else
    {
      if(u > (-b^2/(2*a)))
      {
        val <- Inf
        return(val)
      }
      else
      {
        val <- -b/a - sqrt((b^2+2*a*u)/(a^2))
        return(val)
      }
    }
  }
}

FirstEventTime2 <- function(x)
{
  u <- -log(runif(1))
  val <- log(exp(x) + u/s0) - x
  return(val)
}

FirstEventTime3 <- function(lambda)
{
  return(-log(runif(1))/lambda)
}

EventProb_zz <- function(v,x,t,k)
{
  val1 <- -Y0[k]*v[k] + sum(Sigma_Inv[k,]*(x-mu0))*v[k] + sum(Sigma_Inv[k,]*v)*t*v[k]
  val2 <- s0*exp(x[k]+t*v[k])*v[k]
  
  valTrue <- max(c(0,val1+val2)) + lambda0
  valUpper <- max(c(0,val1)) + max(c(0,val2)) + lambda0
  return(valTrue/valUpper)
}

ZZ <- function(x0,v0)
{
  Niter_zz <- 1e3
  Xtrace_zz <- matrix(0,nrow=Niter_zz, ncol=d)
  Vtrace_zz <- matrix(0,nrow=Niter_zz, ncol=d)
  Ttrace_zz <- rep(0,Niter_zz)
  timeRecord_zz <- matrix(0,nrow=Niter_zz,ncol=d)
  eventRecord_zz <- matrix(0,nrow=Niter_zz,ncol=d)
  
  Xtrace_zz[1,] <- x0
  Vtrace_zz[1,] <- v0
  
  time_zz <- proc.time()
  for(i in 2:Niter_zz)
  {
    tau <- rep(0,d)
    for(k in 1:d)
    {
      a0 = sum(Sigma_Inv[k,]*Vtrace_zz[i-1,]) * Vtrace_zz[i-1,k]
      b0 = Vtrace_zz[i-1,k]*(sum(Sigma_Inv[k,]*(Xtrace_zz[i-1,]-mu0)) - Y0[k])
      time1 <- FirstEventTime1(a0,b0)
      if(Vtrace_zz[i-1,k] == 1)
      {
        time2 <- FirstEventTime2(Xtrace_zz[i-1,k])
      }
      else
      {
        time2 <- Inf
      }
      time3 <- FirstEventTime3(lambda0)
      tau[k] <- min(c(time1, time2, time3))
      timeRecord_zz[i,k] <- tau[k]
      eventRecord_zz[i,k] <- which.min(c(time1, time2, time3))
    }
    Ind <- which.min(tau)
    Ttrace_zz[i] <- tau[Ind]
    
    q <- EventProb_zz(Vtrace_zz[i-1,],Xtrace_zz[i-1,],Ttrace_zz[i],Ind)
    u <- runif(1)
    
    if(u > q)
    {
      Xtrace_zz[i,] <- Xtrace_zz[i-1,] + Ttrace_zz[i]*Vtrace_zz[i-1,]
      Vtrace_zz[i,] <- Vtrace_zz[i-1,]
    }
    else
    {
      Xtrace_zz[i,] <- Xtrace_zz[i-1,] + Ttrace_zz[i]*Vtrace_zz[i-1,]
      Vtrace_zz[i,] <- Vtrace_zz[i-1,]
      Vtrace_zz[i,Ind] <- -Vtrace_zz[i,Ind]
    }
    #print(c(i, Xtrace_zz[i,1:2], Ttrace_zz[i]))
  }
  
  time_zz <- proc.time() - time_zz
  
  fact_zz <- 0.025
  Num_zz <- floor(sum(Ttrace_zz)/fact_zz)
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
  
  return(list(Xsamp=Xsamp_zz, Vsamp=Vsamp_zz, timeHist=timeRecord_zz,
              eventHist = eventRecord_zz, time=time_zz, timeTrace=Ttrace_zz))
}

d <- d0^2
V0 <- sample(c(1,-1),size=d,replace = TRUE)
Result1 <- ZZ(X0,V0)

plot(Result1$Xsamp[,2],type="l")
Result1$time

ResultZZ <- list()
ResultZZ[[1]] <- Result1
for(k in 2:10)
{
  N0 <- nrow(ResultZZ[[k-1]]$Xsamp)
  ResultZZ[[k]] <- ZZ(ResultZZ[[k-1]]$Xsamp[N0,], ResultZZ[[k-1]]$Vsamp[N0,])
  print(c(k,ResultZZ[[k]]$time,N0))
}

XsampZZ <- NULL
for(k in 1:10)
{
  XsampZZ <- rbind(XsampZZ, ResultZZ[[k]]$Xsamp)
  print(c(k, nrow(ResultZZ[[k]]$Xsamp)))
}

save(ResultZZ,file="ResultZZ.RData")

lp_ZZ <- -apply(XsampZZ,1,U)

nrow(XsampCS)
nrow(XsampZZ)

k <- 1
plot(XsampCS[,k],type="o",col="red",pch=20)
lines(XsampZZ[,k],type="o",col="blue",pch=20)

par(mfrow=c(1,2))
par(mar=c(4,4,2,2))
plot(XsampCS,pch=16,ylim=c(-2,6),xlim=c(-2,6),
     xlab=expression(x[1]),ylab=expression(x[2]),main='CS')
plot(XsampZZ,pch=16,ylim=c(-2,6),xlim=c(-2,6),
     xlab=expression(x[1]),ylab=expression(x[2]),main='ZZ')

par(mfrow=c(1,2))
plot(lp_CS,type="l",pch=20,col="red",xlab="",ylab="lp",lwd=2)
lines(lp_ZZ, type="l",pch=20,col="blue",lwd=2)
plot(XsampCS[,400],type="l",col="red",xlab="",ylab=expression(x[400]),lwd=2)
lines(XsampZZ[,400],type="l",col="blue",lwd=2)
par(mfrow=c(1,1))

