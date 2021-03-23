lag_selection<- function(RR, lg) {

  # calculation of lag selection criteria (Akaike, Hannan-Quinnn and Schwarz)
  # the smaller is the criterion, the more adequate is the number of lags
  #
  # INPUT
  # RR: residuals from an estimated VAR
  # lg: number of lags
  # 
  # OUTPUT
  # value for Akaike, Hannan-Quinnn and Schwarz criteria

  T <- nrow(RR)
  K <- ncol(RR)
  M <- RR[1,] %*% t(RR[1,]) # sum of squared errors for first variable
  for (j in 2:T){
    M <- M + (RR[j,] %*% t(RR[j,]))
  }
  M <- M / T # average total SSE

  # AKAIKE criterion
  AK <- log(det(M)) + (2 *lg* K^2)/T
  # Hannan-Quinn criterion
  HQ <- log(det(M)) + (2*log(log(T)) *lg* K^2)/T
  # Schwarz criterion
  SC <- log(det(M)) + (log(T) *lg* K^2)/T
  
  sel <- c(AK, HQ, SC)
  sel
}
