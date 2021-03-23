#' Get names for bilinear effects.
#' @keywords internal
setupMultVarNames  <- function(mult_vars){
  
  if(is.null(mult_vars)){ return(NULL) }
  
  ml <- strsplit(mult_vars, "*", fixed = TRUE)
  
  namesMult <- vector(mode="character", length=length(mult_vars))
  
  for(i in 1:length(namesMult)){namesMult[i] <- paste0(ml[[i]][1],"_by_", ml[[i]][2])}
  
  return(namesMult)
  
}



