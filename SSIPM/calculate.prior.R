#' Calculate the (log) prior given a fitted parameter value
#'
#' @return Value
#'
#' @export
#
# 
calculate.prior <- function(Fit,Prior){
  
  # could include error code if length(Fit) != length(Prior)
  
  Result = rep(NA,length(Fit))
  
  for (i in length(Fit)){
    if (Prior$Type[i] == 'lognormal'){
      
      Result[i] = log(dnorm(log(Fit[i]),mean = Prior$Means[i],sd = Prior$SDs[i]))
  
                      } # end if lognormal (currently there are no other options)
  } # end loop over Fit
  
  Result = sum(Result)
  return(Result)
}
