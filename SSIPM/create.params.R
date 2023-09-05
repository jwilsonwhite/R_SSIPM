#' Generate list of parameter values for the SSIPM
#'
#' @return List
#'
#' @export
#
# 
create.params <- function(Sp,meshsize=100,MCMClen = 100,MCMCchains =2, regulation = 1){
  
  # Other inputs may be necessary in the future
  if (Sp == 'PACL'){ 
    if(regulation == 1){ #account for change in fishing regulations in 2013
    
    # NB I used dummy numbers at first
    M = 0.18 # natural mortality (y^-1)
    Linf = 69.8 # von Bert asymptote, cm
    k = 0.06 # vB growth rate, y^-1
    Lfish = 30.5 # length of 50% recruitment to fishery (cm) - #Catch regulations year from survey start for fishery 
    Lmat = 22.3 # length of 50% maturity (cm)
    Lvar = 0.24 # CV of length (dimensionless)
    Rmean = 5.89 # average size of new recruits
    Rsd = 1.53 # sd of recruit size distribution
    
    # IPM details
    x = seq(0,Lfish*2,length.out=meshsize) # grid of length intervals
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
    if (regulation == 2){
      M = 0.18 # natural mortality (y^-1)
      Linf = 69.8 # von Bert asymptote, cm
      k = 0.06 # vB growth rate, y^-1
      Lfish = 36.5 # length of 50% recruitment to fishery (cm) - #Catch regulations year from survey start for fishery 
      Lmat = 22.3 # length of 50% maturity (cm)
      Lvar = 0.24 # CV of length (dimensionless)
      Rmean = 5.89 # average size of new recruits
      Rsd = 1.53 # sd of recruit size distribution
      
      # IPM details
      x = seq(0,Lfish*2,length.out=meshsize) # grid of length intervals
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
  }
  else if( Sp == "PANE"){
    if(regulation == 1){
      
      M = 0.22 # natural mortality (y^-1)
      Linf = 60.6 # von Bert asymptote, cm
      k = 0.09 # vB growth rate, y^-1
      Lfish = 30.5 # length of 50% recruitment to fishery (cm)
      Lmat = 22.9 # length of 50% maturity (cm)
      Lvar = 0.24 # CV of length (dimensionless)
      Rmean = 5.89 # average size of new recruits
      Rsd = 1.53 # sd of recruit size distribution
      
      # IPM details
      x = seq(0,Lfish*2,length.out=meshsize) # grid of length intervals
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
    if (regulation == 2) {
      M = 0.22 # natural mortality (y^-1)
      Linf = 60.6 # von Bert asymptote, cm
      k = 0.09 # vB growth rate, y^-1
      Lfish = 36.5 # length of 50% recruitment to fishery (cm)
      Lmat = 22.9 # length of 50% maturity (cm)
      Lvar = 0.24 # CV of length (dimensionless)
      Rmean = 5.89 # average size of new recruits
      Rsd = 1.53 # sd of recruit size distribution
      
      # IPM details
      x = seq(0,Lfish*2,length.out=meshsize) # grid of length intervals
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
  }else if (Sp == "SEPU"){
    
    M = 0.2 # natural mortality (y^-1)
    Linf = 56.9 # von Bert asymptote, cm
    k = 0.146 # vB growth rate, y^-1
    Lfish = 36.4 # length of 50% recruitment to fishery (cm)
    Lmat = 24 # length of 50% maturity (cm)
    Lvar = 0.11 # CV of length (dimensionless)
    Rmean = 5.12 # average size of new recruits
    Rsd = 2.14 # sd of recruit size distribution
    
    # IPM details
    x = seq(0,Lfish*2,length.out=meshsize) # grid of length intervals
    dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
    
    # Recruit size distribution
    Rvec = dnorm(x,mean=Rmean,sd=Rsd)
    
    # other model parameters
    burnin = 20 # how many years to initialize the model pre-data
    Q = 20 # now many particles in the filter
    
    # MCMC details
    MCMClen = MCMClen # now many MCMC iterations
    MCMCchains = MCMCchains # how many chains
    
    }else{
      stop('Unrecognized species name')
      }
  
 # Create list for output
  fix.parm <- list(Linf,k,M,Lfish,Lmat,Lvar,Rvec,x,dx,burnin,meshsize,Q,MCMClen,MCMCchains)
  names(fix.parm) <- c('Linf','k','M','Lfish','Lmat','Lvar',
                       'Rvec','x','dx','burnin','meshsize','Q','MCMClen','MCMCchains')
  
 return(fix.parm) 
}
