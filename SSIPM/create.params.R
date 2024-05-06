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
    Fmean = 34.9 # average size of harvest
    Fsd = 4.32 # sd of harvest size 
    Immean = 27.0 # average size of immigrant
    Imsd = 9.96 #sd of immigrant size
    
    # IPM details
    x = seq(0,Linf*2,length.out=meshsize) # grid of length intervals
    dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
    
    # Recruit size distribution
    Rvec = dnorm(x,mean=Rmean,sd=Rsd)
    
    # Fishing selectivity
    F.sel = pnorm(x,mean = Fmean, sd = Fsd)
    
    # Immigrant size distribution
    Ifish = length(x[x<9]) #size of Immigration cutoff (cm)
    Ivec = dlnorm(x, meanlog = log(Immean^2 / sqrt(Imsd^2 + Immean^2)), sdlog = sqrt(log(1 + (Imsd^2 /Immean^2))))
    
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
      Fmean = 36.6 # average size of harvest
      Fsd = 5.67 # sd of harvest size 
      Immean = 27.0 # average size of immigrant
      Imsd = 9.96 #sd of immigrant size
      
      # IPM details
      x = seq(0,Linf*2,length.out=meshsize) # grid of length intervals
      dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
      
      # Recruit size distribution
      Rvec = dnorm(x,mean=Rmean,sd=Rsd)
      
      # Fishing selectivity
      F.sel = pnorm(x,mean = Fmean, sd = Fsd)
      
      # Immigrant size distribution
      Ifish = length(x[x<9]) #size of Immigration cutoff (cm)
      
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
      Fmean = 37 # average size of harvest
      Fsd = 5.41 # sd of harvest size 
      Immean = 30.8 # average size of immigrant
      Imsd = 9.55 #sd of immigrant size
      
      # IPM details
      x = seq(0,Linf*2,length.out=meshsize) # grid of length intervals
      dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
      
      # Recruit size distribution
      Rvec = dnorm(x,mean=Rmean,sd=Rsd)
      
      # Fishing selectivity
      F.sel = pnorm(x,mean = Fmean, sd = Fsd)
      
      # Immigrant size distribution
      Ifish = length(x[x<9]) #size of Immigration cutoff (cm)
      
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
      Fmean = 40 # average size of harvest
      Fsd = 5.6 # sd of harvest size 
      Immean = 30.8 # average size of immigrant
      Imsd = 9.55 #sd of immigrant size
      
      # IPM details
      x = seq(0,Linf*2,length.out=meshsize) # grid of length intervals
      dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
      
      # Recruit size distribution
      Rvec = dnorm(x,mean=Rmean,sd=Rsd)
      
      # Fishing selectivity
      F.sel = pnorm(x,mean = Fmean, sd = Fsd)
      
      # Immigrant size distribution
      Ifish = length(x[x<9]) #size of Immigration cutoff (cm)
      
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
    Fmean = 38.2 # average size of harvest
    Fsd = 8.6 # sd of harvest size 
    Immean = 30.3 # average size of immigrant
    Imsd = 9.96 #sd of immigrant size
    
    # IPM details 
    x = seq(0,Linf*2,length.out=meshsize) # grid of length intervals
    dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
    
    # Recruit size distribution
    Rvec = dnorm(x,mean=Rmean,sd=Rsd)
    
    # Fishing selectivity
    F.sel = pnorm(x,mean = Fmean, sd = Fsd)
    
    # Immigrant size distribution
    Ifish = length(x[x<9]) #size of Immigration cutoff (cm)
    Ivec = dlnorm(x, meanlog = log(Immean^2 / sqrt(Imsd^2 + Immean^2)),
                  sdlog = sqrt(log(1 + (Imsd^2 /Immean^2))))
    
    # other model parameters
    burnin = 20 # how many years to initialize the model pre-data
    Q = 20 # now many particles in the filter
    
    # MCMC details
    MCMClen = MCMClen # now many MCMC iterations
    MCMCchains = MCMCchains # how many chains
    
  }else if (Sp == "CHPU"){
    
    M = 0.65 # natural mortality (y^-1)
    Linf = 22.1 # von Bert asymptote, cm
    k = 0.368 # vB growth rate, y^-1
    Lfish = NA # length of 50% recruitment to fishery (cm)
    Lmat = 14 # length of 50% maturity (cm)
    Lvar = 0.07 # CV of length (dimensionless)
    Rmean = 4.46 # average size of new recruits
    Rsd = 1.51 # sd of recruit size distribution
    Fmean = NA # average size of harvest
    Fsd = NA # sd of harvest size 
    Immean = 15.2 # average size of immigrant
    Imsd = 4.46 #sd of immigrant size
    
    # IPM details 
    x = seq(0,Linf*2,length.out=meshsize) # grid of length intervals
    dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
    
    # Recruit size distribution
    Rvec = dnorm(x,mean=Rmean,sd=Rsd)
    
    # Immigrant size distribution
    Ifish = length(x[x<9]) #size of Immigration cutoff (cm)
    
    #fishing selectivity
    F.sel = NA
    
    # other model parameters
    burnin = 20 # how many years to initialize the model pre-data
    Q = 20 # now many particles in the filter
    
    # MCMC details
    MCMClen = MCMClen # now many MCMC iterations
    MCMCchains = MCMCchains # how many chains
    
    }else if (Sp == "OXCA"){
      
      M = 0.72 # natural mortality (y^-1)
      Linf = 22.45 # von Bert asymptote, cm
      k = 0.564 # vB growth rate, y^-1
      Lfish = NA # length of 50% recruitment to fishery (cm)
      Lmat = 11 # length of 50% maturity (cm)
      Lvar = 0.06 # CV of length (dimensionless)
      Rmean = 4.41 # average size of new recruits
      Rsd = 1.68 # sd of recruit size distribution
      Fmean = NA # average size of harvest
      Fsd = NA # sd of harvest size 
      Immean = 14.0 # average size of immigrant
      Imsd = 4.73 #sd of immigrant size
      
      # IPM details 
      x = seq(0,Linf*2,length.out=meshsize) # grid of length intervals
      dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
      
      # Recruit size distribution
      Rvec = dnorm(x,mean=Rmean,sd=Rsd)
      
      # Immigrant size distribution
      Ifish = length(x[x<9]) #size of Immigration cutoff (cm)
      
      #fishing selectivity
      F.sel = NA
      
      # other model parameters
      burnin = 20 # how many years to initialize the model pre-data
      Q = 20 # now many particles in the filter
      
      # MCMC details
      MCMClen = MCMClen # now many MCMC iterations
      MCMCchains = MCMCchains # how many chains
      
    }else if (Sp == "EMJA"){
      
      M = 0.18 # natural mortality (y^-1)
      Linf = 21.8 # von Bert asymptote, cm
      k = 0.36 # vB growth rate, y^-1
      Lfish = NA # length of 50% recruitment to fishery (cm)
      Lmat = 15 # length of 50% maturity (cm)
      Lvar = 0.13 # CV of length (dimensionless)
      Rmean = 10.77 # average size of new recruits
      Rsd = 2.11 # sd of recruit size distribution
      Fmean = NA # average size of harvest
      Fsd = NA # sd of harvest size 
      Immean = 20.2 # average size of immigrant
      Imsd = 4.73 #sd of immigrant size
      
      # IPM details 
      x = seq(0,Linf*2,length.out=meshsize) # grid of length intervals
      dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
      
      # Recruit size distribution
      Rvec = dnorm(x,mean=Rmean,sd=Rsd)
      
      # Immigrant size distribution
      Ifish = length(x[x<14]) #size of Immigration cutoff (cm)
      
      #fishing selectivity
      F.sel = NA
      
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
  fix.parm <- list(Linf,k,M,Lfish,Lmat,Lvar, Ifish, Rvec, F.sel,x,dx,burnin,meshsize,Q,MCMClen,MCMCchains)
  names(fix.parm) <- c('Linf','k','M', 'Lfish','Lmat','Lvar',
                       'Ifish','Rvec','F.sel', 'x','dx','burnin','meshsize','Q','MCMClen','MCMCchains')
  
 return(fix.parm) 
}
