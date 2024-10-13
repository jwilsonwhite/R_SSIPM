#' DIC model selection of SSIPM
#'
#'
#' @return List
#'
#' @export
#' 
#' DIC has 2 components - Goodness of fit (Deviance) + complexity
#'
#' Complexity  = effective number of parameters (Pd)
#' Deviance = -2*log-likelihood
#' Pd ==‘posterior mean deviance − deviance of posterior means’

ssipm.dic <- function(Data, fix.param, postproc.MCMC){

Comp = list() #list to add components to
  
  # Posterior mean of deviance
  Dbar <- mean(-2*postproc.MCMC$PosteriorProb)

  # deviance of posterior means
  
  Dhat <- var(-2*postproc.MCMC$PosteriorProb)

  #calculate Pd
  Pd <- (Dbar) - (Dhat)


#DIC
DIC <- Dhat + 2*Pd

Comp$Dbar = Dbar
Comp$Dhat = Dhat
Comp$Pd = Pd
Comp$DIC = DIC

return(Comp)

}
