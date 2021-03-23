# PRUNE Prune weight matrix B.
#
# [Bpruned stde ci] = prune(X, k)
#     Prunes and re-estimates the weight matrix using a simple 
#     resampling method. Also re-estimates the standard deviations of error
#     terms and constant terms.
#
#
# Input arguments:
#     X        - the original data matrix used in model estimation.
#     k        - the estimated causal ordering of the variables.
#
# Output:
#     res$Bpruned  - the pruned matrix of connection strenghts.
#     res$stde     - the re-estimated standard deviations of error terms.
#     res$ci       - the re-estimated constant terms.
#
# See also ESTIMATE, LINGAM.

prune <- function( X, k ) {

  # ---------------------------------------------------------------------------
  # Default values for parameters
  # ---------------------------------------------------------------------------

  method <- 'resampling'  # the pruning method

  # 'prunefactor' determines how easily weak connections are pruned
  # in the simple resampling based pruning. Zero gives no pruning 
  # whatsoever. We use a default of 1, but remember that it is quite arbitrary
  # (similar to setting a significance threshold). It should probably
  # become an input parameter along with the data, but this is left
  # to some future version.
  prunefactor <- 1

  # ---------------------------------------------------------------------------
  # Pruning
  # ---------------------------------------------------------------------------

  cat('Pruning the network connections...\n')
  dims <- nrow(X)
  ndata <- ncol(X)

  if (method == 'resampling') {

    # -------------------------------------------------------------------------
    # Pruning based on resampling: divide the data into several equally
    # sized pieces, calculate B using covariance and QR for each piece
    # (using the estimated causal order), and then use these multiple
    # estimates of B to determine the mean and variance of each element.
    # Prune the network using these.
  
    npieces <- 10
    piecesize <- floor(ndata/npieces)

    Bpieces <- array(0,c(dims,dims,npieces))
    diststdpieces <- array(0,c(dims,npieces))
    cpieces <- array(0,c(dims,npieces))
    Bfinal <- matrix(0,dims,dims)
    
    for (i in 1:npieces) {

      # Select subset of data, and permute the variables to the causal order
      Xp <- X[k,((i-1)*piecesize+1):(i*piecesize)]
    
      # Remember to subract out the mean 
      Xpm <- rowMeans(Xp)
      Xp <- Xp - Xpm
    
      # Calculate covariance matrix
      C <- (Xp %*% t(Xp))/ncol(Xp)
    
      # Do QL decomposition on the inverse square root of C
      res <- tridecomp(solve(sqrtm(C)),'ql');
      Q <- res$A
      L <- res$B
    
      # The estimated disturbance-stds are one over the abs of the diag of L
      newestdisturbancestd <- 1/diag(abs(L));
    
      # Normalize rows of L to unit diagonal
      L <- L/diag(L)
    
      # Calculate corresponding B
      Bnewest <- diag(dims)-L;

      # Also calculate constants
      cnewest <- L %*% Xpm

      # Permute back to original variable order
      ik <- as.vector(iperm(k));
      Bnewest <- Bnewest[ik, ik]
      newestdisturbancestd <- newestdisturbancestd[ik]
      cnewest <- cnewest[ik]

      # Save results
      Bpieces[,,i] <- Bnewest
      diststdpieces[,i] <- newestdisturbancestd
      cpieces[,i] <- cnewest

    }

    if (npieces>1) {

      for (i in 1:dims) {
        for (j in 1:dims) {
      
          themean <- mean(Bpieces[i,j,])
          thestd <- sd(Bpieces[i,j,])
          if (abs(themean)<prunefactor*thestd) {
            Bfinal[i,j] <- 0
          }
          else {
            Bfinal[i,j] <- themean
          }
        
        }
      }

    }

    else Bfinal <- Bpieces[,,1] # used when prunefactor = 0
 
    diststdfinal <- rowMeans(diststdpieces)
    cfinal <- rowMeans(cpieces)

    # Finally, rename all the variables to the way we defined them
    # in the function definition
  
    Bpruned <- Bfinal
    stde <- diststdfinal
    ci <- cfinal

  }
  
  if (method == 'olsboot') {
    stop('Not implemented yet!')
  }
  
  if (method == 'wald') {
    stop('Not implemented yet!')
  }

  if (method == 'bonferroni') {
    stop('Not implemented yet!')
  }

  if (method == 'hochberg') {
    stop('Not implemented yet!')
  }

  if (method == 'modelfit') {
    stop('Not implemented yet!')
  }

  cat('Done!\n')

  # Return the result
  res <- list()
  res$Bpruned <- Bpruned
  res$stde <- stde
  res$ci <- ci
  res

}

