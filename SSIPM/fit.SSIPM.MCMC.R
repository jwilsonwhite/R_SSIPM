#' Fit SSIPM to data using delayed rejection MCMC
#'
#' @return List
#'
#' @export
#
# 
fit.SSIPM.MCMC <- function(fix.param,Data,Prior,savename,CVs = c(1,1/2,1/3,1/4,1/5)){
  
  # Set up MCMC
  M = fix.param$MCMClen
  Chains = fix.param$MCMCchains
  
  
  # pre-allocated results variable
  mc_str = list()
  
  for (c in 1:Chains){ # if possible it would be nice to parallelize this step.
    
    ChainP <- rep(0,length.out=M) # this holds the posterior proportional evidence (LL + Prior)
    Values <- matrix(NA,nrow=M,ncol=length(Prior)) 
  
  # Get initial candidate parameter vector  
    cand.param <- get.cand(Values[1,],Prior[c],index=NA,CVs[1])
    Fit <- run.IPM(fix.param,cand.param,Data) # returns list Fit with log-likelihood and fit to data
    Prior.tmp <- calculate.prior(Fit,Prior) # calculate the prior
    ChainP[1] <- Fit$L + Prior.tmp
    
    # parameter counter for one-at-a-time
    k = 1
    # counter for delayed rejection
    kk = 1
    advance = FALSE
    
    for (m in 1:M){ # MCMC steps
      
      while (!advance){  #delayed rejection
        
      # Generate candidate parameters (one-at-a-time)
        #### MAKE TYPE OF CANDIDATE A VARIABLE, OR PUT IT IN THE PRIOR
      cand.param <- get.cand(Values[m,],Prior,type=k,CVs[kk])
      
  # Run the state-space IPM
  Fit <- run.IPM(fix.param,cand.param,Data) # returns list Fit with log-likelihood and fit to data
  Prior.tmp <- calculate.prior(Fit,Prior) # calculate the prior
  Evidence <- Fit$L + Prior.tmp
  
  # Metropolis-Hastings
  p = min(1,exp(Evidence - ChainP[m-1])) # Metropolis step, calculating acceptance probability (note this assumes the candidate generating function is symmetric)
  Accept = runif(1) < p # accept the proposal
  
  if (Accept) { # proposal accepted
    Values[m,] = cand.param
    ChainP[m] = Evidence
    advance = TRUE
    kk = 1
    
  } else if (!Accept & kk < length(Prior)) {  
    kk = kk + 1 # advance to the next value of the CV
  } else {# else if we are out of delayed rejection steps
    # chain stays put
    Values[m,] = Values[m-1,]
    ChainP[m] = ChainP[m-1]
    advance = TRUE
    kk = 1
    
  } # end if Accept
      } # end delayed rejection while loop
      advance = FALSE
  
    } # end loop over m MCMC iterations
    
    mc_str[[c]]$ChainP <- ChainP
    mc_str[[c]]$Values <- Values
    
  } # end loop over chains
  save(mc_str,file=savename)
return(mc_str)
}
