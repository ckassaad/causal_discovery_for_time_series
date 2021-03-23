# randnetbalanced - create a more balanced random network
#
# INPUT:
#
# dims       - number of variables
# indegree   - number of parents of each node (Inf = fully connected)
# parminmax  - [min max] standard deviation owing to parents 
# errminmax  - [min max] standard deviation owing to error variable
#
# OUTPUT:
# 
# res$B      - the strictly lower triangular network matrix
# res$errstd - the vector of error (disturbance) standard deviations
#

randnetbalanced <- function( dims, indegree, parminmax, errminmax ) {

  # Number of samples used to estimate covariance structure
  samples <- 10000 
    
  # First, generate errstd
  errstd <- runif(dims)*(errminmax[2]-errminmax[1]) + errminmax[1]

  # Initializations
  X <- matrix(0,dims,samples)
  B <- matrix(0,dims,dims)

  # Go trough each node in turn:
  for (i in 1:dims) {
   
    # If 'indegree' is finite, randomly pick that many parents,
    # else, all previous variables are parents
    if (!is.infinite(indegree)) {
	    if (i<=indegree) {
          if (i>1) {
	          par <- 1:(i-1)
          } else {
            par <- rep(0,0)
          }
	    } else { 
          par <- sample(1:(i-1),indegree)
      }
    } else {
      if (i>1) {
	      par <- 1:(i-1)
      } else {
        par <- rep(0,0)
      }
    }
    # If node has parents, do the following
    if (length(par)>0) {
	
	# Randomly pick weights
	w <- rnorm(length(par))
	wfull <- matrix(0,i-1,1)
        wfull[par] <- w

	# Calculate contribution of parents
	X[i,] <- t(wfull) %*% X[1:(i-1),]
		
	# Randomly select a 'parents std' 
	parstd <- runif(1)*(parminmax[2]-parminmax[1]) + parminmax[1]
	
	# Scale w so that the combination of parents has 'parstd' std
	scaling <- parstd/sqrt(mean(X[i,]^2));
	w <- w*scaling

	# Recalculate contribution of parents
	wfull <- matrix(0,i-1,1)
        wfull[par] <- w	
	X[i,] <- t(wfull) %*% X[1:(i-1),]
	
	# Fill in B
	B[i,par] <- t(w)

      }
    else {
    # if node has no parents
	
      # Increase errstd to get it to roughly same variance
      parstd <- runif(1)*(parminmax[2]-parminmax[1]) + parminmax[1]
      errstd[i] <- sqrt(errstd[i]^2 + parstd^2)
	
      # Set data matrix to empty
      X[i,] = matrix(0,1,samples)
	
    }
	
    # Update data matrix
    X[i,] = X[i,] + matrix(rnorm(samples),1,samples)*errstd[i]
    
  }

  # Return result
  res <- list()
  res$B <- B
  res$errstd <- errstd
  res
  
}
