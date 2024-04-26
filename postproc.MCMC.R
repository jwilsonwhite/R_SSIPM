#' Calculate the (log) prior given a fitted parameter value
#'
#' @importFrom invgamma dinvgamma
#'
#' @return Value
#'
#' @export
#
# 
postproc.MCMC <- function(Fit){
  
  Chains = length(Fit$ChainP)
  
  Post = list()
  var_within = rep(NA,Chains)
  All.ChainP = c()
  All.Values = c()
  names(All.Values) = names(Fit$Values[[1]])
  
  for (c in 1:Chains){
  # Loop over chains, remove burnin, stack each parameter
  # Calculate values needed for rhat
    
    ChainP = Fit$ChainP[[c]]
    Values = Fit$Values[[c]]
    
    Len = length(ChainP)
    Burnin = round(Len/2)
    
    # remove Burnin
    ChainP = ChainP[Burnin:Len]
    Values = Values[Burnin:Len,]
    
    # Calculate variance
    var_within[c] = var(ChainP)
    
    # Concatenate results
    All.ChainP = c(All.ChainP,ChainP)
    All.Values = rbind(All.Values,Values)
    
  } # end loop over chains
  
  var_comb = var(All.ChainP)
  
  Rhat <- sqrt( var_comb / mean(var_within) ) # Rhat statistic
  
  # Assemble list of output variables
  Post$PosteriorProb = All.ChainP
  Post$Values = All.Values
  Post$Rhat = Rhat
  
  return(Post)
}