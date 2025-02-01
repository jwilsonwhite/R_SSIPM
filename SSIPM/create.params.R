#' Generate list of parameter values for the SSIPM
#'
#' @return List
#'
#' @export
#
# 
create.params <- function(Sp,MCMClen = 100,MCMCchains =2, ogive = NULL, regulation = NULL){
  
  # Other inputs may be necessary in the future
  if (Sp == 'PACL'){ 
    
    # Life history parameters
    M = 0.18 # natural mortality (y^-1)
    Linf = 69.8 # von Bert asymptote, cm
    k = 0.06 # vB growth rate, y^-1
    Lfish = 30.5 # length of 50% recruitment to fishery (cm) - #Catch regulations year from survey start for fishery 
    Lmat = 22.3 # length of 50% maturity (cm)
    Lvar = 0.24 # CV of length (dimensionless)
    Rmean = 5.89 # average size of new recruits
    Rsd = 1.53 # sd of recruit size distribution
    Immean = 27.0 # average size of immigrant
    Imsd = 9.96 #sd of immigrant size
    
    # Fishing before & after regulation change in 2013
    Lfish1 = 30.5 # length of 50% recruitment to fishery (cm) - #Catch regulations year from survey start for fishery 
    Fmean1 = 34.9 # average size of harvest
    Fsd1 = 4.32 # sd of harvest size 
    
    Lfish2 = 36.5 # length of 50% recruitment to fishery (cm) - #Catch regulations year from survey start for fishery 
    Fmean2 = 36.6 # average size of harvest
    Fsd2 = 5.67 # sd of harvest size 

    
    # IPM details
    # x = seq(0,Linf*2,length.out=meshsize) # grid of length intervals
    x = 1:(ceiling(Linf*2))
    meshmax = ceiling(Linf*2)
    dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
    
    # Recruit size distribution
    Rvec = dnorm(x,mean=Rmean,sd=Rsd)
    
    # Fishing selectivity
    F.sel = NULL
    F.sel1 = pnorm(x,mean = Fmean1, sd = Fsd1) #fishing selectivity before 2013
    F.sel2 = pnorm(x,mean = Fmean2, sd = Fsd2) #fishing selectivity 2013 & after
    
    # Immigrant size distribution
    Ifish = length(x[x<9]) #size of Immigration cutoff (cm)
    Ivec = NULL
    
    #minimum observable size
    ogive = 9
    
    # other model parameters
    burnin = 20 # how many years to initialize the model pre-data
    Q = 100 # now many particles in the filter
    
    # MCMC details
    MCMClen = MCMClen # now many MCMC iterations
    MCMCchains = MCMCchains # how many chains

    
  }
  else if( Sp == "PANE"){
      
    # Life history parameters
      M = 0.22 # natural mortality (y^-1)
      Linf = 60.6 # von Bert asymptote, cm
      k = 0.09 # vB growth rate, y^-1
      Lfish = 30.5 # length of 50% recruitment to fishery (cm)
      Lmat = 22.9 # length of 50% maturity (cm)
      Lvar = 0.24 # CV of length (dimensionless)
      Rmean = 5.89 # average size of new recruits
      Rsd = 1.53 # sd of recruit size distribution
      Immean = 30.8 # average size of immigrant
      Imsd = 9.55 #sd of immigrant size
      
      # Fishing before & after regulation change in 2013
      Lfish1 = 30.5 # length of 50% recruitment to fishery (cm)
      Fmean1 = 37 # average size of harvest
      Fsd1 = 5.41 # sd of harvest size 
      
      Lfish2 = 36.5 # length of 50% recruitment to fishery (cm)
      Fmean2 = 40 # average size of harvest
      Fsd2 = 5.6 # sd of harvest size 
      
      # IPM details
      # x = seq(0,Linf*2,length.out=meshsize) # grid of length intervals
      x = 1:(ceiling(Linf*2))
      meshmax = ceiling(Linf*2)
      dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
      
      # Recruit size distribution
      Rvec = dnorm(x,mean=Rmean,sd=Rsd)
      
      # Fishing selectivity
      F.sel = NULL
      F.sel1 = pnorm(x,mean = Fmean1, sd = Fsd1) #fishing selectivity before 2013
      F.sel2 = pnorm(x,mean = Fmean2, sd = Fsd2) #fishing selectivity 2013 & after
      
      # Immigrant size distribution
      Ifish = length(x[x<9]) #size of Immigration cutoff (cm)
      
      #minimum observable fish size
      ogive <- 9
      
      # other model parameters
      burnin = 20 # how many years to initialize the model pre-data
      Q = 100 # now many particles in the filter
      
      # MCMC details
      MCMClen = MCMClen # now many MCMC iterations
      MCMCchains = MCMCchains # how many chains
    
  }else if (Sp == "SEPU"){
    
  # Life history parameters
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
    # x = seq(0,Linf*2,length.out=meshsize) # grid of length intervals
    x = 1:(ceiling(Linf*2))
    meshmax = ceiling(Linf*2)
    dx = diff(x[1:2]) # grid width (needed for midpoint rule integration)
    
    # Recruit size distribution
    Rvec = dnorm(x,mean=Rmean,sd=Rsd)
    # Rvec.tmp[length(x[x<9]):meshsize] = 0 #size cutoff for YOY (recruit)
    # Rvec = Rvec.tmp/sum(Rvec.tmp*fix.param$dx) #standardize probability density
    # 
    # Fishing selectivity
    F.sel = pnorm(x,mean = Fmean, sd = Fsd) #no change in fishing regulation
    F.sel1 = NULL
    F.sel2 = NULL
    
    # Immigrant size distribution
    Ifish = length(x[x<9]) #size of Immigration cutoff (cm)
    Ivec = NULL
    
    #minimum observable fish size
    ogive <- 1
    
    # other model parameters
    burnin = 20 # how many years to initialize the model pre-data
    Q = 100 # now many particles in the filter
    
    # MCMC details
    MCMClen = MCMClen # now many MCMC iterations
    MCMCchains = MCMCchains # how many chains
    
  }else{
      stop('Unrecognized species name')
      }
  
 # Create list for output
  fix.parm <- list(Linf,k,M,Lfish,Lmat,Lvar, Ifish, Rvec, F.sel,F.sel1,F.sel2,ogive, x, meshmax, dx,burnin,meshsize,Q,MCMClen,MCMCchains)
  names(fix.parm) <- c('Linf','k','M', 'Lfish','Lmat','Lvar',
                       'Ifish','Rvec','F.sel','F.sel1','F.sel2', 'ogive', 'x','meshmax', 'dx','burnin','meshsize','Q','MCMClen','MCMCchains')
  
 return(fix.parm) 
}
