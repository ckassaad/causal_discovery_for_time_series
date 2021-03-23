co_ml <- function(Y, p, r){

  # Maximum Likelihood estimation of a cointegrated VAR process (Johansen 
  # procedure) see Lütkepohl, Introduction to Multiple Time series, pp.356-357
  #
  # INPUT
  # Y: original data matrix (not in canonical form)
  # p: number of lags
  # r: cointegration rank
  # 
  # OUTPUT
  # list containing the following
  #   $const: contant term of the process
  #   $residuals: residuals
  #   $D: <-DD
  #   $Pi: <-H%*%C
  #   $AA: VAR level coefficient matrices


  dims <- dim(Y)
  n <- dims[1]
  k <- dims[2]

  # get diferences
  DYt <- Y[2:n,]-Y[1:(n-1),]
  DYs <- t(DYt)
  DY <- DYs[,p:(n-1)]

  T <- n-p
  DXs <- matrix(nrow=(k*(p-1)), ncol=T)
  cont <- 1
  for(i in 1:(p-1)){
    for(j in 1:k){
      DXs[cont,] <- DYs[j,(p-i):(n-1-i)]
      cont <- cont+1
    }
  }
  on <- rep(1,T)
  DX <- rbind(on, DXs)
  Y_pt <- Y[1:(n-p),]
  Y_p <- t(Y_pt)
  I_T <- diag(T)
  M <- I_T-(t(DX)%*%solve(DX%*%t(DX))%*%DX)
  R_0 <- DY%*%M
  R_1 <- Y_p%*%M
  S_00 <- (R_0%*%t(R_0))/T
  S_10 <- (R_1%*%t(R_0))/T
  S_01 <- (R_0%*%t(R_1))/T
  S_11 <- (R_1%*%t(R_1))/T
  P <- t(chol(S_11))
  G <- solve(P)

  L <- G %*% S_10 %*% solve(S_00) %*% S_01 %*% t(G)
  # first r eigenvectors already sorted according decreasing eigenvalues
  ev <- eigen(L)$vectors[,1:r]
  C <- t(ev)%*%G
  H <- -S_01%*%t(C)%*%solve(C%*%S_11%*%t(C))
  D <- (DY+H%*%C%*%Y_p)%*%t(DX)%*%solve(DX%*%t(DX))
  resi <- t(DY-D%*%DX+H%*%C%*%Y_p)
  Sigma_u <- (t(resi)%*%resi)/T

  ## VAR level matrices of coefficients ###
  DD <- D[,-1]
  AA <- list()
  AA[[1]]<- diag(k) +   DD[,1:k]
  for (i in 2:(p-1)) {
    AA[[i]] <- DD[,(i*k-k+1):(i*k)] - DD[,((i-1)*k-k+1):((i-1)*k)]
  }
  AA[[p]] <- - DD[,(i*k-k+1):(i*k)] - (H%*%C)

  # output
  results <- list()
  results$const <- D[,1]
  results$residuals <- resi
  results$D <- DD
  results$Pi <- H%*%C
  results$AA <- AA
  results

}
