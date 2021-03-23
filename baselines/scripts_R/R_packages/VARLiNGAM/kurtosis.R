kurtosis <- function( X ) {

  # INPUT
  # X: data vector or data matrix with variables in columns
  # 
  # OUTPUT
  # kurtosis of the (normalized) columns of X

  X <- as.matrix(X)
  n <- ncol(X)
  kurt <- array(0,dim=c(1,n))

  for (i in 1:n) {

    # normalize X to zero mean and unit variance
    norm <- (X[,i]-mean(X[,i]))/sqrt(var(X[,i]))

    # calculuate excess kurtosis
    kurt[i] <- mean(norm^4)-3*(mean(norm^2))^2

  }

  kurt

}
