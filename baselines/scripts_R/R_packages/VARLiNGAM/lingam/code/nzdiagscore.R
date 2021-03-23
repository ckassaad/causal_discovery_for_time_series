nzdiagscore <- function( W ) {

  res <- sum(1/diag(abs(W)))
  res
  
}
