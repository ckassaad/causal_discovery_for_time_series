MA_components <- function (AA,maxlag){

  # calculate the Moving Average (MA) representation of a VAR model
  #
  # INPUT
  # AA: list of estimated coefficients of reduced VAR
  #     Ahat[[1]] = A_1, ... Ahat[[p]] = A_p
  # maxlag: highest time lag to be included in analysis
  # 
  # OUTPUT
  # a list of matrices containing the MA coefficients up to 'maxlag' time lags

  k <- nrow(AA[[1]]) # number of variables
  p <- length(AA) # time lages

  # define coefficients matrices Ahat_q=0 for p<q<=maxlag
  AA_n <- as.list(1:maxlag)
  for(i in 1:maxlag) {
    AA_n[[i]] <- matrix(0,nrow=k,ncol=k)
  }
  for(i in 1:p) {
    AA_n[[i]] <- AA[[i]]
  }

  # moving average representation
  # use recursion: Fi_0 = I_k
  #                Fi_i = sum_j=1..i Fi_(i-j) AA_n_j, i=1,2,...,maxlag
  Fi <- as.list(1:(maxlag+1))
  for(i in 1:(maxlag+1)) {
    Fi[[i]] <- matrix(0,nrow=k,ncol=k)
  }
  Fi[[1]] <- diag(k) # Fi_0 = I_k
  Fi[[2]] <- AA_n[[1]] # Fi_1 = Fi_0 AA_1 by recursion formula
  
  for (i in 2:maxlag) {
    # recursion formula
    for (j in 1:i) {
      FF <- Fi[[i+1-j]]%*%AA_n[[j]]
      Fi[[i+1]] <- Fi[[i+1]] + FF
    }
  }

  Fi
}


