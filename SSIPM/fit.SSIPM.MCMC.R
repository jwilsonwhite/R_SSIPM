#' Fit SSIPM to data using delayed rejection MCMC
#'
#' @return List
#'
#' @export
#
# 
fit.SSIPM.MCMC <- function(fix.param,Data,Prior,savename,CVs = c(1,1/2,1/3,1/4,1/5),
                           burnin.im = TRUE, burnin.r = TRUE, ipm.im = TRUE, ipm.r = TRUE){
  
  # Set up MCMC
  M = fix.param$MCMClen
  Chains = fix.param$MCMCchains
  
  # pre-allocated results variable
  mc_str = list()
  
  for (c in 1:Chains){ # if possible it would be nice to parallelize this step.
    
    ChainP <- rep(0,length.out=M) # this holds the posterior proportional evidence (LL + Prior)
    Values <- matrix(NA,nrow=M,ncol =length(Prior$Names)) 
    
    Values[1,] <- exp(Prior$Means)
    colnames(Values) <- Names
  
  # Get initial candidate parameter vector  
    cand.param <- get.cand(Values[1,],Prior,index=NA,CVs[1])

    # initiate Immigration size distribution 
    if(isTRUE(ipm.im)|isTRUE(burnin.im)){ #if immigration is used
      im.cand.param <-  c( rep(1,51), 0, cand.param["F"],0) #use updated cand.param F value
      names(im.cand.param) <- c(paste0("r",c(0:50)), "Im",'F','error')
      im.data <- Data
      im.data[,c(1:50)] <- NA #not reliant on data
      im.ipm <- run.IPM(fix.param, im.cand.param, im.data, burnin.im = FALSE, ipm.im = FALSE)
      stbstate <- im.ipm$N[,50] # Distant future based on IPM
      stbstate[c(1:fix.param$Ifish)] <- 0
      
      #replace Ivec with new stable state distribution
      fix.param$Ivec <- stbstate
    } # end if immigration is used
    
    Fit <- run.IPM(fix.param,cand.param,Data, burnin.im = burnin.im, burnin.r = burnin.r, ipm.im = ipm.im, ipm.r = ipm.r) # returns list Fit with log-likelihood and fit to data
    Prior.tmp <- calculate.prior(cand.param,Prior) # calculate the prior
    ChainP[1] <- Fit$LL + Prior.tmp
    Values[1,] <- cand.param 
    
    # parameter counter for one-at-a-time
    k = 1
    # counter for delayed rejection
    kk = 1
    advance = FALSE
    
    for (m in 2:M){ # MCMC steps
      
      while (!advance){  #delayed rejection
        
      # Generate candidate parameters (one-at-a-time)  
      cand.param <- get.cand(Values[m-1,],Prior,index=k,CVs[kk]) 
      
      # update Immigration size distribution at each step
      if(isTRUE(ipm.im )| isTRUE(burnin.im)){ #if immigration is used
      im.cand.param <-  c( rep(1,51), 0, cand.param["F"],0) #use updated cand.param F value
      names(im.cand.param) <- c(paste0("r",c(0:50)), "Im",'F','error')
      im.data <- Data
      im.data[,c(1:50)] <- NA #not reliant on data
      im.ipm <- run.IPM(fix.param, im.cand.param, im.data, burnin.im = FALSE, ipm.im = FALSE)
      stbstate <- im.ipm$N[,50] # Distant future based on IPM
      stbstate[c(1:fix.param$Ifish)] <- 0
      
      #replace Ivec with new stable state distribution
      fix.param$Ivec <- stbstate
      }
      
  # Run the state-space IPM
  Fit <- run.IPM(fix.param,cand.param,Data, burnin.im = burnin.im, burnin.r = burnin.r, ipm.im = ipm.im, ipm.r = ipm.r) # returns list Fit with log-likelihood and fit to data
  Prior.tmp <- calculate.prior(cand.param,Prior) # calculate the prior
  Evidence <- Fit$LL + Prior.tmp
  
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
      # move on to next parameter
      if (k < length(cand.param)){
      k = k + 1}
      else{ k = 1}
      
    } # end loop over m MCMC iterations
    
    mc_str$ChainP[[c]] <- ChainP
    mc_str$Values[[c]] <- Values

  } # end loop over chains
  save(mc_str,file=savename)
return(mc_str)
}
