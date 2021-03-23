nzdiagbruteforce <- function( W ) {

  #--------------------------------------------------------------------------
  # Try all row permutations, find best solution
  #--------------------------------------------------------------------------

  n <- nrow(W)
  
  bestval <- Inf;
  besti <- 0;
  allperms <- all.perm(n) 
  nperms <- nrow(allperms)
  
  for (i in 1:nperms) {
    Pr <- diag(n)
    Pr <- Pr[,allperms[i,]]
    Wtilde <- Pr %*% W
    c <- nzdiagscore(Wtilde)
    if (c<bestval) {
      bestWtilde <- Wtilde
      bestval <- c
      besti <- i
    }
  }

  res <- list()
  res$Wopt <- bestWtilde
  res$rowp <- allperms[besti,]
  res$rowp <- iperm(res$rowp)
  
  res
  
}
