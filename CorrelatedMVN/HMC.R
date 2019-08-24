library(MASS)
library(mcmcse)
library(parallel)

HMC_Sampler <- function(epsilon, Niter, d, seed0)
{
  set.seed(seed0 * 1e7)
  
  A <- matrix(0,nrow=d,ncol=d)
  for(i in 1:d)
  {
    for(j in 1:d)
      A[i,j] <- 0.9^(abs(i-j))
  }
  A_Inv <- solve(A)
  
  U <- function(x)
  {
    return(as.numeric(0.5 * t(x)%*% A_Inv %*% x))
  }
  
  grad_U <- function(x)
  {
    return(as.vector(A_Inv %*% x))
  }
  
  Leapfrog <- function(U,grad_U,epsilon,L,current_q)
  {
    ### Based on the Neal's implementation ###
    q = current_q
    p = rnorm(length(q))
    current_p = p
    p = p - 0.5*epsilon * as.vector(grad_U(q))
    
    for(l in 1:L)
    {
      q = q + epsilon * as.vector(p)
      if(l != L)
        p = p - epsilon * as.vector(grad_U(q))
    }
    p = p - 0.5*epsilon * grad_U(q)
    p = -p
    current_U <- U(current_q)
    current_K <- 0.5*as.numeric(sum(current_p^2))
    proposed_U <- U(q)
    proposed_K <- 0.5*as.numeric(sum(p^2))
    rho <- exp(current_U-proposed_U+current_K-proposed_K)
    if(is.na(rho))
    {
      return(c(current_q,0,NA))
    }
    else
    {
      if(runif(1) < rho)
      {
        return(c(q,1,rho))
      }
      else
      {
        return(c(current_q,0,rho))
      }
    }
  }
  
  Niter_hmc <- Niter
  Xsamp_hmc <- matrix(0,nrow=Niter_hmc,ncol=d)
  Xsamp_hmc[1,] <- mvrnorm(n=1,mu=rep(0,d),Sigma=A)
  Comp <- 0
  ratio <- rep(0,Niter_hmc)
  ratio_theory <- rep(0,Niter_hmc)
  
  # step size for leapfrog
  epsilon <- epsilon
  time_hmc <- proc.time()
  for(i in 2:Niter_hmc)
  {
    # number of steps for leapfrog integrator
    r <- sample(1:45,size=1)
    Xstar <- Leapfrog(U,grad_U,epsilon,r,Xsamp_hmc[(i-1),])
    Xsamp_hmc[i,] <- Xstar[1:d]
    ratio[i] <- Xstar[d+1]
    Comp <- Comp + r
    if(!is.na(Xstar[d+2]))
    {
      ratio_theory[i] <- Xstar[d+2]
    }
    else
    {
      ratio_theory[i] <- 0
    }
    #if(i %% floor(Niter_HMC/10) == 0)
    #  print(paste("Iteration ",toString(100*round(i/Niter_HMC,digits = 2)), "%",sep=""))
  }
  time_hmc <- proc.time() - time_hmc
  #print(c(sum(ratio)/Niter_HMC,time_HMC))
  
  ks_hmc <- rep(0,d)
  for(i in 1:d)
  {
    ks_hmc[i] <- as.numeric(ks.test(jitter(Xsamp_hmc[,i],amount=0.00001), "pnorm", 0, 1)[[1]])
  }
  
  return(list(Xsamp_ess=ess(Xsamp_hmc), time=time_hmc, grad_num=Comp,
              Xsamp_ks = ks_hmc))
  
}

d0 <- 20
Niter0 <- 4e5

num_cores <- detectCores()

cl <- makeCluster(num_cores,type = "FORK")

ResultHMC <- parLapply(cl, 1:40, function(zz) HMC_Sampler(0.3,Niter0,d0,zz))

save(ResultHMC, file="ResultHMC.RData")
stopCluster(cl)

