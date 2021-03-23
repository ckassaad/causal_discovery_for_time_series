#' Recode variable names.
#' @param data The vector of variable names to be recoded
#' @param oldvalue A vector containing the latent variable names used internally.
#' @param newvalue A vector containing the observed variable names, either
#' provided by the user (as a header) or provided by R (e.g., V1, V2).
#' @return Recoded vector of variable names.
#' @keywords internal 
recode.vars <- function(data,
                        oldvalue,
                        newvalue){
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  newvec <- data
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  newvec
}
