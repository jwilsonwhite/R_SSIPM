#' Calculate the (log) prior given a fitted parameter value
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
      
      Result[i] = log(dnorm(log(Cand[i]),mean = Prior$Means[i],sd = Prior$SDs[i]))
  
                      } # end if lognormal (currently there are no other options)
  } # end loop over Cand
  
  Result = sum(Result)
  return(Result)
}
