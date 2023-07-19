#' Generate candidate parameter values
#'
#' @return Value
#'
#' @export
#
# 
get.cand <- function(Init,Prior,index = NA,CV = 1){
 
  if (!is.na(index)){ # for just one of the list of parameters
    
    if (Prior.Type[index] == 'lognormal'){
      
      cand = exp( rnorm(n = 1, mean=log(Init), sd = Prior$SDs[index]*CV) )
      # Note that this generates a random variable on the original scale, not log-transformed
      
    } # end if lognormal (currently the only option)
    
  } # end if !is.na
  else{ # then loop over all the values
    
    cand  <- rep(NA,length(Init))
    for (k in 1:length(Init)){
      
      if (Prior.Type[k] == 'lognormal'){
        
        cand[k] = exp( rnorm(n = 1, mean=log(Init[k]), sd = Prior$SDs[k]*CV) )
        # Note that this generates a random variable on the original scale, not log-transformed
        
      } # end if lognormal (currently the only option)
      
    } # end loop over all
    
  } # end else
  
return(cand)   
}