#' Counts number of excellent fit indices
#' @param indices A vector of fit indices from lavaan.
#' @return The number of fit indices that are excellent.
#' @keywords internal 
count.excellent <- function(indices){
  rmseaE    <- ifelse(indices[4] < .05, 1, 0)
  srmrE     <- ifelse(indices[5] < .05, 1, 0)
  nnfiE     <- ifelse(indices[6] > .95, 1, 0)
  cfiE      <- ifelse(indices[7] > .95, 1, 0)
  excellent <- sum(rmseaE, srmrE, nnfiE, cfiE, na.rm = TRUE)
  return(excellent)
}