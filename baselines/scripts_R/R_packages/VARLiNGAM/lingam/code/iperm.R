iperm <- function( p ) {

  q <- array(0,c(1,length(p)))
  
  for (i in 1:length(p)) {
    ind <- which(p==i)
    q[i] <- ind[1]
  }
  
  q
  
}
