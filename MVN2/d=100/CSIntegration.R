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


CS_sampler <- function(TimeGap, Niter, d, seed0)
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
  A_Inv <- solve(A) # inverse of A
  
  T_delta <- TimeGap
  Niter_cs <- Niter
  
  # store the sample points of (x,v)
  Xsamp_cs <- matrix(0,nrow=Niter_cs,ncol=d)
  Vsamp_cs <- matrix(0,nrow=Niter_cs,ncol=d)
  
  # initialize (x_1, v_1)
  Xsamp_cs[1,] <- mvrnorm(n=1,mu=rep(0,d),Sigma=A)
  Ind <- sample(1:d,size=1)
  Vsamp_cs[1,Ind] <- sample(c(-1,1),size=1)

  # record the total virtual time along the path
  T_sum <- 0
  
  # remark the current (x,v)
  X_current <- Xsamp_cs[1,]
  V_current <- Vsamp_cs[1,]
  A_Inv_x <- as.vector(A_Inv %*% X_current)
  i <- 0
  
  # grad_num : number of computing gradient
  # Lambda_num: number of recall FirstEventTime
  grad_num <- 1 
  Lambda_num <- 0
  
  time_cs <- proc.time()
  while(i < Niter_cs)
  {
    T_Increment <- FirstEventTime(A_Inv[Ind,Ind], V_current[Ind]*A_Inv_x[Ind])
    # increment Lambda_num by 1
    Lambda_num <- Lambda_num + 1
    X_propose <- X_current + T_Increment * V_current
    A_Inv_x <- as.vector(A_Inv %*% X_propose)
    # increment grad_num
    grad_num <- grad_num + 1
    p = abs(A_Inv_x)
    p = p/sum(p)
    Ind <- sample(1:d,size=1,prob=p)
    V_propose <- rep(0,d)
    if(A_Inv_x[Ind] > 0)
      V_propose[Ind] = -1
    else
      V_propose[Ind] = 1
    T_sum_new <- T_sum + T_Increment
    Inte_delta <- floor(T_sum_new/T_delta) - i
    if(Inte_delta >= 1)
    {
      for(n in 1:Inte_delta)
      {
        i <- i + 1
        if(i <= Niter_cs)
        {
          Xsamp_cs[i,] <- X_current + (i*T_delta - T_sum)  * V_current
          Vsamp_cs[i,] <- V_current
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
  time_cs <- proc.time() - time_cs
  return(list(Xsamp_ess = ess(Xsamp_cs), time=time_cs, Lambda_num=Lambda_num, grad_num=grad_num))
}

TimeGap0 <- 0.1
d0 <- 100
Niter0 <- 5e3

num_cores <- detectCores()

cl <- makeCluster(num_cores,type = "FORK")

ResultCS <- parLapply(cl, 1:40, function(zz) CS_sampler(d0,Niter0,d0,zz))

save(ResultCS, file="ResultCS.RData")
stopCluster(cl)


