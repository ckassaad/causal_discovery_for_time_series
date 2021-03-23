irf <- function(AA, B0, lag_irf, u_res) {

  # calculate impulse response functions
  # 
  # INPUT
  # AA: list of estimated coefficients of reduced VAR
  #     Ahat[[1]] = A_1, ... Ahat[[p]] = A_p
  # B0: instantaneous effects
  # lag_irf: number of lags included in the calculation of the irf
  # u_res: residuals from reduced form VAR
  # 
  # OUTPUT
  # IRF: impulse response functions in k^2 x lag_irf matrix, with k=#variables
  #      each row contains the response of one variable on another varialbe at
  #      all time lags (0,...,lag_irf)
  #      rows 1...k contain the response of the 1st var to all other vars
  #      rows k+1...2k contain the response of the 2nd var to all other vars
  #      etc.

  k <- ncol(u_res) # number of variables

  # structural residuals
  Gamma0 <- diag(k)-B0 # = P^(-1)
  v_res <- t(Gamma0 %*% t(u_res))

  # Moving average (MA) representation of the reduced VAR process
  FI <- MA_components(AA,lag_irf)
  P <- solve( Gamma0 ) # Gamma0 = P^(-1)
  PSI <- FI
  for (i in 1:length(FI)) {
    # Moving average representation of the structural VAR = Fi * P
    PSI[[i]] <- FI[[i]] %*% P
  }

  # impulse response functions
  IRF <- matrix(nrow=k*k, ncol=lag_irf+1)
  count <- 1
  for (i in 1:k) {
    for (j in 1:k) {
      for (g in 1:(lag_irf+1)) {
        # irf of lag g, response of variable i to a unit shock in variable j
        IRF[count,g] <- PSI[[g]][i,j]*sd(v_res[,j])
      }
    count<-count+1
    }
  }

  IRF
}
