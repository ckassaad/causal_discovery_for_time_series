sltscore <- function( B ) {

  s <- sum((B[upper.tri(B,diag=TRUE)])^2)
  s
  
}
