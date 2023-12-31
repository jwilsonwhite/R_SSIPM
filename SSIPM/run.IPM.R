#' Run a state-space IPM
#'
#' @return List
#'
#' @export
#
# 
run.IPM <- function(fix.param,cand.param,Data){
  
  params = c(fix.param,cand.param)
  
  # Create the kernel
  K = kernmat(params,timestep=1) #
  
  # Initialize the model:
  datayears = dim(Data)[2] # how many years of observations
  
  N0 = matrix(0,nrow=params$meshsize,ncol=params$burnin) # get the initial distribution
  N0[,1] = params$Rvec * params$r0 # initialize with one pulse of new recruits
  
  # Run the model for the burnin period
  for (t in 2:params$burnin){
    N0[,t] <- K %*% N0[,t-1] * params$dx  + params$r0 * params$Rvec # midpoint rule integration
  }
  
  # Now run for the time period when we have data
  N = matrix(0,nrow=params$meshsize,ncol=datayears+1)
  N[,1] = N0[,params$burnin]
  L = rep(NA,datayears)
  
  for (t in 2:datayears){
    
    #advance the model
    Rt = get("params") [[paste0("r",t)]] # extract the recruitment for this year
    N[,t] = K %*% N[,t-1] * params$dx  + Rt * params$Rvec
  
    # apply particle filter & calculate likelihood
    if (!is.na(Data[1,t])){ # if there are observations in this model year, otherwise we just skip it
    Nq = matrix(0,nrow=params$meshsize,ncol=params$Q)
    Lt = rep(NA,params$Q)
    
    for (q in 1:params$Q){
  
      rlm <- N[,t]
      rls <- params$error
      Nq[,q] = rlnorm(n=params$meshsize,meanlog=log(rlm^2 / sqrt(rls^2 + rlm^2))
                      ,sdlog = sqrt(log(1 + (rls^2 /rlm^2)))) # lognormal process error
  
      Lt[q] = sum(dpois( x = Data[[t]], lambda = Nq[,q], log = TRUE))# poisson likelihood (note - requires integer count data
      # Could include a correction here for the number of transects observed, if that differs
    }
    
    # reweight the observations based on likelihood
    # advance the weighted average
    mLt = sum(Lt)
    Ltt = Lt/mLt
    Ltm = t(matrix(Ltt,nrow=params$Q,ncol=params$meshsize))
    
    Ntmp = rowSums(Nq * Ltm) # weighted average
    
    N[,t] = Ntmp
    
    L[t-1] =  sum( dpois( x = Data[[t]], lambda = N[,t], log = TRUE) ) # the likelihood
    
    } # end if data
    
  } # end loop over data years
  
  Ltotal = sum(L,na.rm= TRUE)
  
  IPM.result = list(LL = Ltotal,N = N)
  
  
 return(IPM.result) 
}
