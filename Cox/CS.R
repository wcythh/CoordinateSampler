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
EventProb <- function(v,x,t,k)
{
  val1 <- -Y0[k]*v + v* sum(Sigma_Inv[k,]*(x-mu0)) + t*Sigma_Inv[k,k]
  val2 <- s0*exp(x[k])*exp(t*v)*v
  
  valTrue <- max(c(0, val1+val2)) + lambda0
  valUpper <- max(c(0, val1)) + max(c(0, val2)) + lambda0
  
  return(valTrue/valUpper)
}

FirstEventTime1 <- function(a,b)
{
  u = -log(runif(1))
  if(b > 0)
    return(sqrt((b^2+2*u*a)/a^2)-b/a)
  else
    return(sqrt(2*u/a)+abs(b/a))
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

lambda0 <- 0.1

CS <- function(x0,v0)
{
  Niter_cs <- 1e4
  Xtrace_cs <- matrix(0,nrow=Niter_cs,ncol=d)
  Vtrace_cs <- matrix(0,nrow=Niter_cs,ncol=d)
  Ttrace_cs <- rep(0,Niter_cs)
  
  Xtrace_cs[1,] <- x0
  Vtrace_cs[1,] <- v0
  
  Ind <- which.max(abs(Vtrace_cs[1,]))
  
  time_cs <- proc.time()
  timeRecord <- matrix(0,nrow=Niter_cs,ncol=3)
  for(i in 2:Niter_cs)
  {
    a0 <- Sigma_Inv[Ind, Ind]
    b0 <- Vtrace_cs[i-1,Ind] *(sum(Sigma_Inv[Ind,]*(Xtrace_cs[i-1,]-mu0)) - Y0[Ind])
    time1 <- FirstEventTime1(a0,b0)
    if(Vtrace_cs[i-1,Ind] == 1)
    {
      time2 <- FirstEventTime2(Xtrace_cs[i-1,Ind])
      time3 <- FirstEventTime3(lambda0)
    }
    else
    {
      time2 <- Inf
      time3 <- FirstEventTime3(lambda0)
    }
    Ttrace_cs[i] <- min(c(time1, time2, time3))
    timeRecord[i,] <- c(time1,time2,time3)
    q <- EventProb(Vtrace_cs[i-1,Ind], Xtrace_cs[i-1,], Ttrace_cs[i], Ind)
    
    u <- runif(1)
    if(u >q)
    {
      Xtrace_cs[i,] <- Xtrace_cs[i-1,] + Ttrace_cs[i] * Vtrace_cs[i-1,]
      Vtrace_cs[i,] <- Vtrace_cs[i-1,] 
    }
    else
    {
      Xtrace_cs[i,] <- Xtrace_cs[(i-1),] + Ttrace_cs[i] * Vtrace_cs[(i-1),]
      gradx <- grad_U(Xtrace_cs[i,])
      
      p1 <- cbind(-gradx, 0)
      p2 <- cbind(gradx,0)
      
      p <- c(apply(p1,1,max), apply(p2,1,max)) + lambda0
      Ind <- sample(1:(2*d), size=1, prob=p)
      if(Ind <= d)
      {
        Vtrace_cs[i,Ind] <- 1
      }
      else
      {
        Ind <- Ind - d
        Vtrace_cs[i,Ind] <- -1
      }
    }
    #print(c(i,time1,time2,time3))
  }
  time_cs <- proc.time() - time_cs
  
  fact_cs <- 100
  Num_cs <- floor(sum(Ttrace_cs)/fact_cs)
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
  return(list(Xsamp=Xsamp_cs, Vsamp=Vsamp_cs, timeHist=timeRecord, time=time_cs))
}

#lp <- -apply(Xsamp_cs,1,U)
#plot(lp,type="l",pch=20,lwd=2)
d <- d0^2
V0 <- rep(0,d)
Ind <- sample(1:d,size=1)
V0[Ind] <- sample(c(1,-1),size=1)

Result1 <- CS(X0,V0)

plot(Result1$Xsamp[,1],type="l")
Result1$time

ResultCS <- list()
ResultCS[[1]] <- Result1
for(k in 2:10)
{
  N0 <- nrow(ResultCS[[k-1]]$Xsamp)
  ResultCS[[k]] <- CS(ResultCS[[k-1]]$Xsamp[N0,], ResultCS[[k-1]]$Vsamp[N0,])
  print(c(k,ResultCS[[k]]$time,N0))
}

XsampCS <- NULL
for(k in 1:10)
{
  XsampCS <- rbind(XsampCS, ResultCS[[k]]$Xsamp)
  print(c(k, nrow(ResultCS[[k]]$Xsamp)))
}

save(ResultCS,file="ResultCS.RData")

lp_CS <- -apply(XsampCS,1,U)

