boot_sd_firmgrowth <- function(Growth, nboot, Ahat, Bhat) {

  # Bootstrap the firm growth data to get standard errors of coefficients of
  # reduced and structural VAR
  # 
  # INPUT
  # Growth: Firmgrowth data in canonical form
  # cons: estimated constant term of reduced VAR
  # Ahat: list of estimated coefficients of reduced VAR
  #       Ahat[[1]] = A_1, ... Ahat[[p]] = A_p
  # Bhat: list of estimated coefficients of structural VAR
  #       Bhat[[1]] = B_0, ... Bhat[[p+1]] = B_p
  # u_res: residuals from reduced VAR
  # ord: causal order (f.ex. found by LiNGAM)
  # p: time lags
  # nboot: number of bootstrap samples
  # 
  # OUTPUT
  # prints coefficients, standard errors, p-values and an indicator of
  # significance at level 0.01 for matrices A_1, ... A_p, and B_0, B_1, ... B_p

  nvar <- 4
  nlags <- (dim(Growth)[2] - 2)/nvar - 1

  dims1 <- dim(Growth)

  MM <- list()
  BB <- list()
  for (i in 1:nlags) {
    MM[[i]] <- array(0,dim=c(nboot,nvar^2))
    BB[[i]] <- array(0,dim=c(nboot,nvar^2))
  }
  BB0 <- array(0,dim=c(nboot,nvar^2))

  start <- proc.time()
  for (i in 1:nboot) {

    cat("bootstrap run", i, "out of", nboot, "\n")
    
    # get bootstrap sample
    pick <- sample(dims1[1],dims1[1],replace=TRUE)
    Growth_boot <- Growth[pick,]
    res <- VARLiNGAM(Growth_boot,"lad", ntests=FALSE)
    
    # write results columnwise in rows
    BB0[i,] <- res$Bhat[[1]]
    for (j in 1:nlags) {
      MM[[j]][i,] <- res$Mhat[[j]]
      BB[[j]][i,] <- res$Bhat[[j+1]]
    }

  }
  end <- proc.time()
  end-start

  # only use those models which have same instantaneous effects than original
  # model
  bool <- BB0[,2]==0 & BB0[,5]!=0 & BB0[,12]!=0 & BB0[,15]==0

  cat("Number of models with same causal order as orignial one: ", sum(bool),
      "\n")

  sdMM <- list()
  sdBB <- list()
  for (i in 1:nlags) {
    sdMM[[i]] <- sd(MM[[i]][bool,])
    dim(sdMM[[i]]) <- c(4,4)
    sdBB[[i]] <- sd(BB[[i]][bool,])
    dim(sdBB[[i]]) <- c(4,4)
  }
  sdBB0 <- sd(BB0[bool,])
  dim(sdBB0) <- c(4,4)

  # get t-statistics
  alpha <- 0.01
  n <- dim(Growth)[1]+nlags
  df <- 4^2*nlags + 4*3/2 + 4 # number of parameters:
    # nlags 4x4 matrices of lagged effects (full matrices estimated)
    # one 4x4 matrix of instantaneous effects (only lower triangle estimated)
    # one 4x1 vector of constant terms
  tvalueBs <- qt(1-alpha/2, n-df-1)
  df <- 4^2*nlags + 4 # number of parameters:
    # nlags 4x4 matrices of lagged effects (full matrices estimated)
    # one 4x1 vector of constant terms
  tvalueMs <- qt(1-alpha/2, n-df-1)


  tstatsMM <- list()
  tstatsBB <- list()
  for (i in 1:nlags) {
    tstatsMM[[i]] <- Ahat[[i]]/sdMM[[i]]
    tstatsBB[[i]] <- Bhat[[i+1]]/sdBB[[i]]
    if (i==1 | i==2) {
      cat('\ninformation for A', i, ': coeffs, sd, pvalue, significant\n')
      print(round(Ahat[[i]],4))
      print(round(sdMM[[i]],4))
      print(2*pt(abs(tstatsMM[[i]]),n,lower.tail=FALSE))
      print(abs(tstatsMM[[i]]) > tvalueMs)
      cat('\ninformation for B', i, ': coeffs, sd, pvalue, significant\n')
      print(round(Bhat[[i+1]],4))
      print(round(sdBB[[i]],4))
      print(2*pt(abs(tstatsBB[[i]]),n,lower.tail=FALSE))
      print(abs(tstatsBB[[i]]) > tvalueBs)
    }
  }
  tstatsBB0 <- Bhat[[1]]/sdBB0
  cat('\ninformation for B0', ': coeffs, sd, pvalue, significant\n')
  print(round(Bhat[[1]],4))
  print(round(sdBB0,4))
  print(2*pt(abs(tstatsBB0),n,lower.tail=FALSE))
  print(abs(tstatsBB0) > tvalueBs)

}