#rm(list=ls(all=TRUE))
library(MASS)
library(mcmcse)
library(parallel)


FirstEventTime <- function(a, b)
{
  ## generate time interval between two events for the rate function with the form {a*t + b}_+
  u = -log(runif(1))
  if(a > 0)
  {
    if(b > 0)
      return((sqrt((b^2+2*a*u)) - b)/a)
    else
      return(sqrt(2*u/a) - b/a)
  }
  if(a == 0)
  {
    if(b > 0)
      return(u/b)
    else
      return(Inf)
  }
  if(a < 0 )
  {
    if(b <= 0)
      return(Inf)
    else
    {
      if(u > -b^2/(2*a))
        return(Inf)
      else
        return(-sqrt((2*a*u+b^2)/(a^2)) - b/a)
    }
  }
}


BPS_sampler <- function(TimeGap, Niter, d, seed0)
{
  ##@parameters:
  ## -- TimeGap: length of vitural time interval between two adjacent sample points
  ## -- Niter:   number of sample points
  ## -- d:       dimension of the target distribution
  ## -- seed0:   random seed for reproduction
  
  # set the random seed for reproduction
  set.seed(seed0 * 1e7)
  
  # set the covariance matrix for multivariate normal distribution
  A <- matrix(0.0,d,d)
  for(i in 1:d)
  {
    for(j in 1:d)
    {
      A[i,j] <- 0.9^(abs(i-j))
    }
  }
  #A <- diag(1, d)
  A_Inv <- solve(A) # inverse of A
  
  T_delta <- TimeGap
  Niter_bps <- Niter
  
  # store the sample points of (x,v)
  Xsamp_bps <- matrix(0,nrow=Niter_bps,ncol=d)
  Vsamp_bps <- matrix(0,nrow=Niter_bps,ncol=d)
  
  # initialize (x_1, v_1)
  Xsamp_bps[1,] <- mvrnorm(n=1,mu=rep(0,d),Sigma=A)
  Vsamp_bps[1,] <- rnorm(d)
  Vsamp_bps[1,] <- Vsamp_bps[1,]/sqrt(sum(Vsamp_bps[1,]^2))
  
  # record the total virtual time along the path
  T_sum <- 0
  
  # remark the current (x,v)
  X_current <- Xsamp_bps[1,]
  V_current <- Vsamp_bps[1,]
  A_Inv_x <- as.vector(A_Inv %*% X_current)
  A_Inv_v <- as.vector(A_Inv %*% V_current)
  i <- 0
  
  # grad_num : number of computing gradient
  # Lambda_num: number of recall FirstEventTime
  grad_num <- 1 
  Lambda_num <- 0
  
  time_bps <- proc.time()
  while(i < Niter_bps)
  {
    A_Inv_v <- as.vector(A_Inv %*% V_current)
    T_Increment <- FirstEventTime(sum(V_current*A_Inv_v), sum(V_current*A_Inv_x))
    # increment Lambda_num by 1
    Lambda_num <- Lambda_num + 1
    X_propose <- X_current + T_Increment * V_current
    A_Inv_x <- as.vector(A_Inv %*% X_propose)
    # increment grad_num
    grad_num <- grad_num + 1
    gradx <- A_Inv_x/sqrt(sum(A_Inv_x^2))
    V1 <- sum(gradx * V_current) * gradx 
    V2 <- V_current - V1
    V_propose <- -V1 + V2
    T_sum_new <- T_sum + T_Increment
    Inte_delta <- floor(T_sum_new/T_delta) - i
    if(Inte_delta >= 1)
    {
      for(nn in 1:Inte_delta)
      {
        i <- i + 1
        if(i <= Niter_bps)
        {
          Xsamp_bps[i,] <- X_current + (i*T_delta - T_sum)  * V_current
          Vsamp_bps[i,] <- V_current
        }
        #print(i)
      }
      X_current <- X_propose
      V_current <- V_propose
      T_sum <- T_sum_new
    }
    else
    {
      X_current <- X_propose
      V_current <- V_propose
      T_sum <- T_sum_new
    }
  }
  time_bps <- proc.time() - time_bps
  
  ks_bps <- rep(0,d)
  for(i in 1:d)
  {
    ks_bps[i] <- as.numeric(ks.test(jitter(Xsamp_bps[,i],amount=0.00001), "pnorm", 0, 1)[[1]])
  }
  
  return(list(Xsamp_ess = ess(Xsamp_bps), time=time_bps, Lambda_num=Lambda_num, grad_num=grad_num, 
              Xsamp_ks = ks_bps))
}

d0 <- 20
Niter0 <- 1.7e5

num_cores <- detectCores()

cl <- makeCluster(num_cores,type = "FORK")

ResultBPS <- parLapply(cl, 1:40, function(zz) BPS_sampler(20,Niter0,d0,zz))

save(ResultBPS, file="ResultBPS.RData")
stopCluster(cl)


