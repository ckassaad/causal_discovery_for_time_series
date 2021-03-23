sltbruteforce <- function( B ) {

  #--------------------------------------------------------------------------
  # Try all row permutations, find best solution
  #--------------------------------------------------------------------------

  n <- nrow(B)
  
  bestval <- Inf;
  besti <- 0;
  allperms <- all.perm(n) 
  nperms <- nrow(allperms)
  
  for (i in 1:nperms) {
    Btilde <- B[allperms[i,],allperms[i,]]
    c <- sltscore(Btilde)
    if (c<bestval) {
      bestBtilde <- Btilde
      bestval <- c
      besti <- i
    }
  }

  res <- list()
  res$Bopt <- bestBtilde
  res$optperm <- allperms[besti,]
  
  res

}
