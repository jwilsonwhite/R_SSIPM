#' Generate list of parameter values for the SSIPM
#'
#' @return List
#'
#' @export
#
# 
create.params <- function(Sp,gridsize=100,MCMClen = 100,MCMCchains =2){
  
  # Other inputs may be necessary in the future
  if (Sp == 'PCLA'){
    
    # NB I used dummy numbers at first
    M = 0.2 # natural mortality (y^-1)
    Linf = 100 # von Bert asymptote, cm
    k = 0.2 # vB growth rate, y^-1
    Lfish = 40 # length of 50% recruitment to fishery (cm)
    Lmat = 40 # length of 50% maturity (cm)
    Lvar = 0.1 # CV of length (dimensionless)
    Rmean = 10 # average size of new recruits
    Rsd = 1 # sd of recruit size distribution
    
    # IPM details
    x = seq(0,Lfish*2,length.out=gridsize) # grid of length intervals
    dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
    
    # Recruit size distribution
    Rvec = dnorm(x,mean=Rmean,sd=Rsd)
    
    # other model parameters
    burnin = 20 # how many years to initialize the model pre-data
    Q = 20 # now many particles in the filter
    
    # MCMC details
    MCMClen = MCMClen # now many MCMC iterations
    MCMCchains = MCMCchains # how many chains
  }
  else{
    stop('Unrecognized species name')
  }
  
  
 # Create list for output
  fix.parm <- list(Linf,k,M,Lfish,Lmat,Lvar,Rvec,x,dx,burnin,Q)
  names(fix.parm) <- c('Linf','k','M','Lfish','Lmat','Lvar',
                       'Rvec','x','dx','burnin','Q')
  
 return(fix.parm) 
}