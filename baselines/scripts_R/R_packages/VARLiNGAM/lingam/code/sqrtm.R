sqrtm <- function( A ) {

  e <- eigen(A)
  V <- e$vectors
  B <- V %*% diag(sqrt(e$values)) %*% t(V)
  B
  
}
