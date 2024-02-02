#' Fit create IPM kernel
#'
#' @return List
#'
#' @export
#
# # Function to create the IPM kernel, given input parameters
kernmat <- function(params,timestep){
  
  # create a mesh grid of size changes
  x = params$x
  X = t(matrix(x,nrow=length(x),ncol=length(x)))
  Y = t(X)
  # X is matrix of sizes at t
  # Y is matrix of sizes at t+1
  
  # Survival part of kernel
  m = params$M
  
  if(all(is.na(params$F.sel))){ #if fishing selectivity is not available
  pm = exp(-m*timestep)
  } else {
    pm = exp(-(m + params$F*params$F.sel )*timestep)
  }
  
  # Growth part of kernel. Do it this way so that we simulate many different growth trajectories in the population, averaged together.
  nLinfs = 1e3 # how many different values of Linf to simulate
  Linfs =rnorm(n = nLinfs, mean=params$Linf,sd=params$Lvar) # vector of distribution of Linfs
  Linfs.mat = matrix(Linfs,nrow=nLinfs,ncol=length(x)) # expand into a matrix so there is a corresponding value of Linf for each possible value in the length vector x
  X.mat = t(matrix(x,ncol=nLinfs,nrow=length(x))) # expand x into a matrix so there is a value for each value of Linfs.mat
  
  g.tmp = Linfs.mat - (Linfs.mat - X.mat)*exp(-params$k*timestep) # use those two matrices to get the range of possible growth rates, as a function of X
  g.mean = colMeans(g.tmp) # Take the mean across all of the different trajectories for each value of x
  g.mat = t(matrix(g.mean,nrow=length(x),ncol=length(x))) # expand into a matrix with a corresponding value for each value of Y (the size at time t+1)
  pg = dnorm(Y,mean=g.mat,sd=2*params$dx) # use dnorm to get the distribution of growth rates (using an arbitrarily small sd)
  
  # make sure no negatives
  pg[pg<0] = 0
  m[m<0] = 0 
  
  # you would add fecundity here if it were a closed population...
  k = pm*pg
  
  # If adding fecundity too:
  #Rvec = dnorm(Y,params$Rmean,params$Rsd); # size distribution of recruits
  #Fec = params$Fec_constant*X^Params$Fec_exponent; # you can use values of 1 and 3 for these... 
  #Q = Fec*Rvec; # combine the two
  #k = k + Q; # add in fecundity to the growth/survival part
  
  
  return(k)
}
