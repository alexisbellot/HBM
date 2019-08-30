## Synthetic data generation for toy example

# This produces survival competing risks data as described in the main paper

simul <- function(N)
{
  require('flexsurv')
  # covariate --> k uniform
  X = runif(N,-2,2)
  
  #latent event times
  Ltime1 = X^3+15 + rnorm(N,0,1) #  1 cubic
  Ltime2 = 10 + 5*X^2 + rnorm(N,0,1) # 2 parabola
  
  
  # censoring times and status
  time = numeric(N)
  ind_cause1 = which(X> -1&X<1)
  u = runif(N,0,1)
  time[ind_cause1] = ifelse(u[ind_cause1]<0.8,Ltime1[ind_cause1],Ltime2[ind_cause1])
  time[-ind_cause1] = ifelse(u[-ind_cause1]<0.8,Ltime2[-ind_cause1],Ltime1[-ind_cause1])
  cause = numeric(N)
  cause[ind_cause1] = ifelse(u[ind_cause1]<0.8,1,2)
  cause[-ind_cause1] = ifelse(u[-ind_cause1]<0.8,2,1)
  
  # data set
  data.frame(id=1:N,
             Survival=time,Ltime1=Ltime1,Ltime2=Ltime2,Cause=cause,
             X=X)
}