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


ZZ_sampler <- function(TimeGap, Niter, d, seed0)
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
  A_Inv <- solve(A)   # inverse of A 
  
  T_delta <- TimeGap
  Niter_zz <- Niter
  
  # store the sample points of (x,v)
  Xsamp_zz <- matrix(0,nrow=Niter_zz,ncol=d)
  Vsamp_zz <- matrix(0,nrow=Niter_zz,ncol=d)
  
  # initialize (x_1, v_1)
  Xsamp_zz[1,] <- mvrnorm(n=1,mu=rep(0,d),Sigma=A)
  Vsamp_zz[1,] <- sample(c(-1,1),size=d,replace = TRUE)
  
  # record the total virtual time along the path
  T_sum <- 0
  
  # remark the current (x,v)
  X_current <- Xsamp_zz[1,]
  V_current <- Vsamp_zz[1,]
  
  X_propose <- X_current
  V_propose <- V_current
  
  # grad_num : number of computing gradient
  # Lambda_num: number of recall FirstEventTime
  grad_num <- 0
  Lambda_num <- 0
  
  i <- 0
  
  time_zz <- proc.time()
  while(i < Niter_zz)
  {
    A_Inv_x <- as.vector(A_Inv %*% X_current)
    A_Inv_v <- as.vector(A_Inv %*% V_current)
    # increment grad_num
    grad_num <- grad_num + 1
    Ttemp <- rep(0,d)
    for(j in 1:d)
      Ttemp[j] <- FirstEventTime(V_current[j]*A_Inv_v[j], V_current[j]*A_Inv_x[j])
    # increment Lambda_num by 1
    Lambda_num <- Lambda_num + d
    Ind <- which.min(Ttemp)
    T_sum_new <- T_sum + Ttemp[Ind]
    X_propose <- X_current + V_current * Ttemp[Ind]
    V_propose <- V_current
    V_propose[Ind] <- - V_propose[Ind]
    Inte_delta <- floor(T_sum_new/T_delta) - i
    
    if(Inte_delta >= 1)
    {
      for(n in 1:Inte_delta)
      {
        i <- i + 1
        if(i <= Niter_zz)
        {
          Xsamp_zz[i,] <- X_current + (i*T_delta - T_sum)  * V_current
          Vsamp_zz[i,] <- V_current
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
  time_zz <- proc.time() - time_zz
  return(list(Xsamp_ess=ess(Xsamp_zz), time=time_zz, Lambda_num=Lambda_num, grad_num=grad_num))
}


TimeGap0 <- 0.1
d0 <- 40
Niter0 <- 5e3

num_cores <- detectCores()

cl <- makeCluster(num_cores,type = "FORK")

ResultZZ <- parLapply(cl, 1:40, function(zz) ZZ_sampler(1,Niter0,d0,zz))

save(ResultZZ, file="ResultZZ.RData")
stopCluster(cl)
