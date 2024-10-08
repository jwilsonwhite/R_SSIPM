#' Run a state-space IPM
#'
#' @return List
#'
#' @export
#
# 
run.IPM <- function(fix.param,cand.param,Data, burnin = TRUE,
                    burnin.im = TRUE, burnin.r = TRUE, ipm.im = TRUE, ipm.r = TRUE){
  
  params = c(fix.param,cand.param)
  
  # Create the kernel
  K = kernmat(params,timestep=1) #
  
  # Initialize the model:
  datayears = dim(Data)[2] # how many years of observations
  
  if(burnin == TRUE){
    N0 = matrix(0,nrow=params$meshmax,ncol=params$burnin) # get the initial distribution
    N0[,1] = params$Rvec * params$r0 # initialize with one pulse of new recruits
    
    # Run the model for the burnin period
    for (t in 2:params$burnin){
      if(isTRUE(burnin.im) & isTRUE( burnin.r)){ #if both r & im during burnin
        N0[,t] <- K %*% N0[,t-1] * params$dx  + params$r0 * params$Rvec + params$Im0 * params$Ivec # midpoint rule integration
      }else
        if(isTRUE(burnin.im) & !isTRUE(burnin.r)){ #if no recruit during burnin
          N0[,t] <- K %*% N0[,t-1] * params$dx + params$Im0*params$Ivec
        }else
          if(isTRUE(burnin.r & !isTRUE(burnin.im))){ #if no recruit during burnin
            N0[,t] <- K %*% N0[,t-1] * params$dx + params$r0 * params$Rvec
          }else
            if(!isTRUE(burnin.im ) & !isTRUE( burnin.r)){
              N0[,t] <- K %*% N0[,t-1] * params$dx 
            }
      
      N0[ N0[,t] <= 1e-323, t ] = 1e-323 #Prevent the generation of Inf values
      
    } #end burnin 
  } #end if run burning

  # Now run for the time period when we have data
  N = matrix(0,nrow=params$meshmax,ncol=datayears) ####!!!!!!!!!1 removed +1 from data years
  
  if(burnin == TRUE){ # if using a burnin period
    N[,1] = N0[,params$burnin]} else { #initialize the model
      if(isTRUE(ipm.im) & isTRUE(ipm.r)){ # if include both im and r term
# browser()
      N[,1] = params$Rvec * params$r1 + params$Ivec*params$Im1
      }
      if(isTRUE(ipm.im) & !isTRUE(ipm.r)){ # if include only im term
      N[,1] =  params$Im1 * params$Ivec
      }
      if(isTRUE(ipm.r) & !isTRUE(ipm.im)){ # if include only r term
      N[,1] =  params$r1 * params$Rvec
      }
      if(!isTRUE(ipm.im) & !isTRUE(ipm.r)) { # if include none
      N[,1] =  0
      }
    } # end if burnin 
  
  colnames(N) = c(1:(datayears)) ####!!!!!!!!!1 removed +1 from data years
  L = rep(NA,datayears)
  N[ N[,1] <= 1e-323, 1 ] = 1e-323 #Prevent the generation of Inf values
  
  for (t in 2:datayears){
    #advance the model
    Rt = get("params") [[paste0("r",t)]] # extract the recruitment for this year
    Im = get("params") [[paste0("Im",params$phase[t])]] # extract immigration for this period
    # Im = get("params") [[paste0("Im",t)]] #!!!!!!!!!!!! TESTING individual immigration
    
    if(isTRUE(ipm.im) & isTRUE(ipm.r)){ # if include both im and r term
      # browser()
      N[,t] = K %*% N[,t-1] * params$dx  + Rt * params$Rvec + Im * params$Ivec 
    }else
      if(isTRUE(ipm.im) & !isTRUE(ipm.r)){ # if include only im term
        N[,t] = K %*% N[,t-1] * params$dx  + Im * params$Ivec
      }else
        if(isTRUE(ipm.r) & !isTRUE(ipm.im)){ # if include only r term
          N[,t] = K %*% N[,t-1] * params$dx  + Rt * params$Rvec
        }else
          if(!isTRUE(ipm.im) & !isTRUE(ipm.r)) { # if include none
            N[,t] = K %*% N[,t-1] * params$dx
          }

    N[ N[,t] <= 1e-323, t ] = 1e-323 #Prevent the generation of Inf values
    
    # apply particle filter & calculate likelihood
    if (!is.na(Data[1,t])){ # if there are observations in this model year, otherwise we just skip it
      Nq = matrix(0,nrow=params$meshmax,ncol=params$Q)
      Lt = rep(NA,params$Q)
      
      for (q in 1:params$Q){
        ###!!!!!!!!! Temporary removed due to issues with calculating process error, use fixed error
        # rlm <- N[,t] # arithmetic expected value
        # rls <- params$error # arithmetic variance
        # rlm[ rlm <= 1e-323] = 1e-323 #Prevent the generation of Inf values
        # meanlog <- log(rlm^2 / sqrt(rls^2 +rlm^2)) # calculate lognormal expected vale
        # if(any(is.infinite(meanlog))){
        # meanlog[meanlog == -Inf] = min(meanlog[meanlog != -Inf])
        # }#end if infinite meanlog
        # sdlog <- sqrt(log(1 + (rls^2 /rlm^2))) #calculate lognormal sd
        # if(any(is.infinite(sdlog))){
        # sdlog[ sdlog == Inf] = max(sdlog[sdlog != Inf])
        # } #end if infinite sdlog
        # Nq[,q] = rlnorm(n=params$meshmax, meanlog=meanlog, sdlog = sdlog) # lognormal process error
        
        Nq[,q] = exp(rnorm(n=params$meshmax, mean= log(N[,t]), sd = params$error)) # lognormal process error
  
        
        Nq[ Nq[,q] <= 1e-323,  ] = 1e-323 #Prevent the generation of Inf values
          Nq.temp = Nq[,q] * params$correction[[1,t]] # correction for differing survey area 
          ####!!!!!!!!!!!!!!!!!!! Convert numerical density to count for posisson estimation
        Data.tmp <- Data[[t]]  ####!!!!!!!!!!!!!!!!!!! Separate data with likelihood calculation
        
        if(length(Data.tmp) != length(Nq.temp)){ #error code if size bins don't match
          stop("Error: Data size bin does not match mesh size bin")
        }
        
        if(!is.null(fix.param$obsize)){ #specify likely observable size
        Lt[q] = sum(dpois( x = Data.tmp[fix.param$obsize:length(Data.tmp)], lambda = Nq.temp[fix.param$obsize:length(Data.tmp)], log = TRUE)) 
        # poisson likelihood (note - requires integer count data)
        } else{
          Lt[q] = sum(dpois( x = Data.tmp, lambda = Nq.temp, log = TRUE)) # poisson likelihood (note - requires integer count data)
          warning("Minimial observation size not specified, likelihood calculated across all sizes")
        }
      }  

      # reweight the observations based on likelihood
      # advance the weighted average
      mLt = sum(Lt)
      Ltt = Lt/mLt
      Ltm = t(matrix(Ltt,nrow=params$Q,ncol=params$meshmax))
      
      Ntmp = rowSums(Nq * Ltm) # weighted average
      
      N[,t] = Ntmp
      
      L[t-1] =  sum( dpois( x = Data[[t]], lambda = N[,t], log = TRUE) ) # the likelihood
  
    } # end if data
  } # end loop over data years
  
  Ltotal = sum(L,na.rm= TRUE)
  
  IPM.result = list(LL = Ltotal,N = N)
  
  return(IPM.result) 
}

