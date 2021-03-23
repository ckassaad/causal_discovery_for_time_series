choleski <- function(u_res, ord) {

  # Choleski factorization on residuals from a reduced form VAR
  # 
  # INPUT
  # u_res: residuals of the reduced form VAR
  # ord: a causal order over the variables
  # nlags: number of time lags
  # corank: cointegration rank
  # 
  # OUTPUT
  # B0: instantaneous effects, in same order as u_res

  u_res_n <- u_res[,ord] # re-ordering the residuals to causal order
  k <- ncol(u_res)

  OM <- cov(u_res_n) # covariance matrix of reduced form residuals 
  P <- t(chol(OM)) # Cholesky factorization of Covariance matrix

  D <- matrix(0,k,k) # scaling matrix to get 1's on diag of W
  diag(D) <- diag(P)

  W <- D %*% solve(P) # rotation matrix
  B0n <- (diag(k)-W) # instantaneous effects
  B0n[upper.tri(B0n, diag=TRUE)] <- 0 # rounding effects
 
  # put back to ordering of the original data
  B0 <- B0n[iperm(ord),iperm(ord)] 
  B0
}

