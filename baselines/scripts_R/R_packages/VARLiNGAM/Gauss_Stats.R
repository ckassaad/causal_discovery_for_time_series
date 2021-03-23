Gauss_Stats <- function(X,n=NULL) {

  # INPUT
  # X: data vector or data matrix (assume more observations than variables)
  # n: number of bins in histograms
  #
  # OUTPUT
  # histograms, qqplot and excess kurtosis of each variable in X

  X <- as.matrix(X)
  dims <- dim(X)

  if (dims[1]<dims[2]) {
    X <- t(X)
    dims <- dim(X)
  }

  if (is.null(n)) {
    if (dims[1] <= 500) n <- 15
    if (dims[1] > 500 & dims[1] <= 2000) n <- 25
    if (dims[1] > 2000) n <- 50
  }

  wi <- dims[2]*2.5
  x11(width=wi, height=5)
  par(mfcol=c(2,dims[2]))
  cat('The excess kurtosis are: \n')

  for (i in 1:dims[2]) {
    # histograms with normal distribution to compare
    hist(X[,i],breaks=n,freq=FALSE,mfg=c(1,i,2,dims[2]))
    m <- mean(X[,i])
    s <- sqrt(var(X[,i]))
    xval <- seq(m-3*s,m+3*s,by=6*s/100)
    # corresponding normal distribution
    lines(xval,dnorm(xval,m,s),col='red',lwd=1.5)

    # qq-plots
    qqnorm(X[,i],mfg=c(2,i,2,dims[2]))

    # normalize data to zero mean and unit variance
    norm <- (X[,i]-m)/s
    # excess kurtosis
    cat('   ', mean(norm^4)-3*(mean(norm^2))^2, '\n')
  }

}