# simpletestlingam - LiNGAM test using completely random parameters
#
# SYNTAX:
# simpletestlingam()
#
# What?
# Randomly selects testing parameter values, generates data, and
# runs LiNGAM to estimate the parameters. Reports and plots the
# results.
#

simpletestlingam <- function() {

  #----------------------------------------------------------------------
  # 1. GENERATE A MODEL (RANDOMLY SELECT PARAMETERS)
  #----------------------------------------------------------------------

  # Number of variables to use: min 2, max 8 (max set arbitrarily)
  dims <- sample(2:8,1)  
  cat('Number of dimensions: ',dims,'\n',sep='')

  # Full or sparse connections?
  sparse <- sample(0:1,1)
  if (sparse == 1) {
    indegree <- sample(1:3,1)    # how many parents maximally 
    cat('Sparse network, max parents = ',indegree,'\n',sep='')
  }
  else { 
    indegree <- Inf              # full connections
    cat('Full network.\n')
  }
  
  parminmax <- c(0.5,1.5) # [min max] standard deviation owing to parents
  errminmax <- c(0.5,1.5) # [min max] standard deviation owing to disturbance

  # Pause briefly to allow user to inspect chosen parameters
  wait(2)

  # Create the network with random weights but according to chosen parameters
  cat('Creating network...')
  res <- randnetbalanced( dims, indegree, parminmax, errminmax )
  B <- res$B
  disturbancestd <- res$errstd
  c <- 2*rnorm(dims) 
  cat('Done!\n')


  #----------------------------------------------------------------------
  # 2. GENERATE DATA FROM THE MODEL
  #----------------------------------------------------------------------

  cat('Generating data...');

  # Nonlinearity exponent, selected to lie in [0.5, 0.8] or [1.2, 2.0].
  # (<1 gives subgaussian, >1 gives supergaussian)
  q <- runif(dims)*1.1+0.5;    
  ind <- which(q>0.8)           
  q[ind] <- q[ind]+0.4     

  # Number of data vectors
  samples <- 10000  

  # This generates the disturbance variables, which are mutually 
  # independent, and non-gaussian
  S <- matrix(rnorm(dims*samples),dims,samples)
  S <- sign(S)*(abs(S)^q);
  
  # This normalizes the disturbance variables to have the 
  # appropriate scales
  S <- S/(sqrt(colMeans(t(S)^2))/disturbancestd)

  # Now we generate the data one component at a time
  Xorig <- matrix(0,dims,samples)
  for (i in 1:dims) {
    Xorig[i,] <- B[i,]%*%Xorig + S[i,] + c[i]
  }
  
  # Select a random permutation because we do not assume that we know
  # the correct ordering of the variables
  p <- sample(1:dims)
  
  # Permute the rows of the data matrix, to give us the observed data
  X <- Xorig[p,]

  # Permute the rows and columns of the original generating matrix B 
  # so that they correspond to the actual data
  Bp <- B[p,p]

  # Permute the generating disturbance stds so that they correspond to
  # the actual data
  disturbancestdp <- disturbancestd[p]

  # Permute the generating constants so that they correspond to
  # the actual data
  cp <- c[p];

  cat('Done!\n')

  #----------------------------------------------------------------------
  # 3. CALL LINGAM TO DO THE ESTIMATION
  #----------------------------------------------------------------------

  res <- lingam(X)
  Best <- res$B
  stde <- res$stde
  ci <- res$ci
  causalperm <- res$k
  
  #----------------------------------------------------------------------
  # 4. PLOT THE ESTIMATED GRAPH AGAINST THE TRUE GRAPH
  #----------------------------------------------------------------------

  # For small dimensions, also display the actual connection matrices
  if (dims<=8) {
    print(Bp)
    print(Best)
  }
  
  # Plotting the graphs
  plotgraph(Bp)
  plotgraph(Best)
  

}
