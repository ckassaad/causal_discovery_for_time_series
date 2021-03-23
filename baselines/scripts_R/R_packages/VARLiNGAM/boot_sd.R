boot_sd <- function(Data, cons, Ahat, Bhat, u_res, ord, p, nboot) {

  # Bootstrap to get standard errors of coefficients of reduced and structural
  # VAR using VECM with cointegration rank 3, and Cholesky
  # 
  # INPUT
  # Data: original Data 
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

  dims <- dim(u_res)
  t <- dims[1] # number of samples
  k <- dims[2] # number of variables

  kurt <- array(0,dim=c(nboot,6)) # to save kurtosis of new residuals
  c <- array(0,dim=c(nboot,6))

  # to save bootstrap results
  MM <- list()
  BB <- list()
  for (i in 1:p) {
    MM[[i]] <- array(0,dim=c(nboot,k^2))
    BB[[i]] <- array(0,dim=c(nboot,k^2))
  }
  BB0 <- array(0,dim=c(nboot,k^2))

  # to get exactly the same results as in paper
  set.seed(1)

  start <- proc.time()
  for (i in 1:nboot) {

    if (i%%10==0) {
      cat("bootstrap run", i, "out of", nboot, "\n")
    }

    # generating the artificial data
    ind <- sample(t, t, replace=TRUE)
    unew <- u_res[ind,] # randomly sample residuals with replacement
    Ynew <- matrix(0,nrow=(t+p),ncol=k)
    Ynew[1:p,] <- as.matrix(Data[1:p,]) # initial time points
    for(ii in (p+1):(t+p)) {
      for(j in 1:p) {
	Y <- Ahat[[j]] %*% t(Data[ii-j,]) # generate new Data
	Ynew[ii,] <- Ynew[ii,] + Y 
      }
      Ynew[ii,] <- cons + Ynew[ii,] + unew[ii-p,]
    }

    Ynew <- as.data.frame(Ynew)
    Data_can_new <- tsdata2canonicalform(Ynew, p)

    # estimate reduced form VAR using a vecm
    res <- VAR_estim(Data_can_new, "vecm", regstats=FALSE, corank=3)

    kurt[i,] <- kurtosis(res$resid) # kurtosis of reduced form VAR
    c[i,] <- res$const

    # write results columnwise in rows of MM[[j]]
    for (j in 1:p) {
      MM[[j]][i,] <- res$Mhat[[j]]
    }

    # calculate instantaneous effects matrix using Cholesky
    B0ch <- choleski(res$resid,ord) # using causal order of original data

    # write results columnwise in rows of BB[[j]]
    BB0[i,] <- B0ch
    Gamma0 <- diag(k) - B0ch
    for (j in 1:p) {
      BB[[j]][i,] <- Gamma0 %*% res$Mhat[[j]]
    }


  }
  end <- proc.time()
  print(end-start)


  # standard errors of bootstrap samples

  sdMM <- list()
  sdBB <- list()
  for (i in 1:p) {
    sdMM[[i]] <- sd(MM[[i]])
    dim(sdMM[[i]]) <- c(k,k)
    sdBB[[i]] <- sd(BB[[i]])
    dim(sdBB[[i]]) <- c(k,k)
  }
  sdBB0 <- sd(BB0)
  dim(sdBB0) <- c(k,k)


  # get t-statistics
  alpha <- 0.01
  n <- dim(Data)[1]
  df <- k^2*p + k*(k-1)/2 + k # number of parameters:
    # p kxk matrices of lagged effects (full matrices estimated)
    # one kxk matrix of instantaneous effects (only lower triangle estimated)
    # one kx1 vector of constant terms
  tvalueBs <- qt(1-alpha/2, n-df-1)
  df <- k^2*p + k # number of parameters:
    # p kxk matrices of lagged effects (full matrices estimated)
    # one kx1 vector of constant terms
  tvalueMs <- qt(1-alpha/2, n-df-1)

  tstatsMM <- list()
  tstatsBB <- list()
  for (i in 1:p) {
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

