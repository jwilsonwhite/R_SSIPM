#' Calculate the (log) prior given a fitted parameter value
#'
#' @importFrom invgamma dinvgamma
#'
#' @return Value
#'
#' @export
#
# 
calculate.prior <- function(Cand,Prior){
  
  # could include error code if length(Cand) != length(Prior)
  if(length(Cand) != length(Prior$Names)) stop("Error: input vectors do not have the same length") #error code
  
  Result = rep(NA,length(Cand))
  
  for (i in 1:length(Cand)){
    if (Prior$Type[i] == 'lognormal'){
      
      Result[i] = log(dnorm(log(Cand[i]) ,mean = Prior$Means[i],sd = Prior$SDs[i]))
  
                      } # end if lognormal 
    
    
    if (Prior$Type[i] == 'invgamma'){
      
      Result[i] = invgamma::dinvgamma(Cand[i], shape = Prior$Means[i], scale = Prior$SDs[i], log = TRUE) 
        # shape and scale same as Prior values for simplicity of code
                      } #end if inverse gamma
    
  } # end loop over Cand
  
  Result = sum(Result)
  return(Result)
}
