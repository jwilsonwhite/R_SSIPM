#' Generate candidate parameter values
#' 
#' @importFrom invgamma rinvgamma
#'
#' @return Value
#'
#' @export



get.cand <- function(Init,Prior,index = NA,CV = 1){

  if (!is.na(index)){ # for just one of the list of parameters
    
    if (Prior$Type[index] == 'lognormal'){
      
      cand = exp( rnorm(n = 1, mean=log(Init)[index], sd = Prior$SDs[index]*CV) )
      # Note that this generates a random variable on the original scale, not log-transformed
      
      Cand.out = Init
      Cand.out[index] = cand
      
    } # end if lognormal 
    
    if (Prior$Type[index] == "invgamma"){
      
      shape =  Prior$Means[index] #constant uninformative prior
      scale = Init[index]*(shape-1)
      
      cand = invgamma::rinvgamma(n = 1, shape = shape, scale = scale) 
      
      Cand.out = Init
      Cand.out[index] = cand
    }#end if inverse gamma
    
  } # end if !is.na
  else{ # then loop over all the values
    
    cand  <- rep(NA,length(Init))
    for (k in 1:length(Init)){
      
      if (Prior$Type[k] == 'lognormal'){
        
        cand[k] = exp( rnorm(n = 1, mean=log(Init[k]), sd = Prior$SDs[k]*CV) )
        # Note that this generates a random variable on the original scale, not log-transformed
        
      } # end if lognormal (currently the only option)
      
      if (Prior$Type[k] == "invgamma"){
        
        shape =  Prior$Means[k] #constant uninformative prior
        scale = Init[k]*(shape-1)
        
        cand[k] = invgamma::rinvgamma(n = 1, shape = shape, scale = scale) 
        
      } #end if inverse gamma
      
    } # end loop over all
    Cand.out = cand
    names(Cand.out) = Prior$Names
  } # end else
return(Cand.out)
}
