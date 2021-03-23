VARLiNGAM <- function(Data, est_meth="ols", ntests=TRUE, pruning=TRUE,
             regstats=FALSE, corank=NA, fmlad=FALSE) {

  # Estimate a VAR-LiNGAM, see
  # - A. Hyvärinen, S. Shimizu, P.O. Hoyer ((ICML-2008). Causal modelling
  #   combining instantaneous and lagged effects: an identifiable model based
  #   on non-Gaussianity; Section 3
  # - A. Hyvärinen, K. Zhang, S. Shimizu, P.O. Hoyer (JMLR-2010). Estimation of
  #   a Structural Vector Autoregression Model Using Non-Gaussianity;
  #   Subsection 5.1.1
  # 
  # -------------------- Step 1: VAR estimation ----------------------------- #
  # -------------------- Step 2: Residuals ---------------------------------- #
  # -------------------- Step 3: LiNGAM on Residuals ------------------------ #
  # ---------------------Step 4: Calculation of Bhats ----------------------- #
  #
  # INPUT
  # Data: time series data in canonical form:
  #       (see function 'tsdata2canonicalform.R')
  # est_meth: string giving the estimation method for the VAR model:
  #       "ols" for for ordinary least squares
  #       "lad" for least absolute deviation (requires package 'quantreg')
  #       "vecm" for an Vector Error Correction Model
  # ntests: boolean, if TRUE tests for normality are carried out, and qq-plots
  #       and histograms are shown (requires packages 'nortest' and 'tseries')
  # pruning: boolean, if TRUE resulting coefficients of LiNGAM (step 3) for
  #       instantaneous effects are pruned using a resampling approach
  # regstats: boolean, only if est_meth = "ols" or "lad" 
  #       if TRUE return statistics for regression coefficients in
  #       VAR-estimation (such as p-values, t-statistics and standard errors)
  # corank: int, only considered when est_meth="vecm", gives the cointegraion
  #        rank in the VECM
  # fmlad: boolean, only considered when est_meth="lad"
  #       if TRUE the fully modified LAD is used (Phillips ET 1995)
  # 
  # OUTPUT
  # list containing:
  #   $resid: residuals of VAR estimation
  #   $const: estimated constant term of VAR
  #   $Mhat: estimated regression coefficients of VAR
  #   $Bhat: estimated regression coefficients of SVAR (adjusted Mhat with
  #          instantenous effects taken account for)
  #   $var_order: a causal order of the variables


  # -------------------------------- step 0 --------------------------------- #

  # number of observed variables
  nvar <- dim(Data)[2] - which(colnames(Data)=="curval1") + 1
  # number of time lags included
  nlags <- (dim(Data)[2] - 2)/nvar - 1

  # some information on the input
  cat("\n------------------------ some information ------------------------\n")
  cat("using",nlags,"time lag(s) and estimating VAR using",est_meth,"method\n")


  # -------------------------------- step 1 --------------------------------- #

  cat("estimate VAR ... ")

  VARres <- VAR_estim(Data, est_meth, regstats, corank, fmlad)

  cat("Done! \n")


  # -------------------------------- step 2 --------------------------------- #

  cat("calculating residuals ... ")

  nhat <- t(VARres$residuals) 
  dims <- dim(nhat)

  cat("Done! \n")


  # -------------------------------- step 3 --------------------------------- #

  if (ntests) {

    cat("\n Histograms, qq-plot and excess kurtosis of the residuals \n ")
    Gauss_Stats(nhat)

    cat("Tests for normality of residuals ... p-values: \n")
    print(Gauss_Tests(nhat))

  }

  cat("\nPerform LiNGAM analysis on residuals ... \n")

  if (pruning) {
    # prune results immediately
    reslg <- lingam(nhat)
    B0hat <- reslg$Bpruned
  }
  else {
    # get order of variables with LiNGAM, calculate matrix B0 with cholesky
    reslg <- estimate(nhat) # LiNGAM

    # get instantaneous effects B0hat using Cholesky facorization
    # arrange residuals in causal order and get covariance matrix 
    OM <- cov(t(nhat)[,reslg$k]) # covariance matrix of reduced form residuals 
    P <- t(chol(OM)) # choleski triangular matrix
    D <- matrix(0,nvar,nvar) # scaling matrix to get 1's on diag of W
    diag(D) <- diag(P)
    W <- D %*% solve(P) # rotation matrix
    B0n <- (diag(nvar) - W)
    B0n[upper.tri(B0n, diag=TRUE)] <- 0 # rounding effects
    # permute back to original variable order
    B0hat <- B0n[iperm(reslg$k),iperm(reslg$k)]
  }

  cat("Done with LiNGAM analysis! \n")


  # -------------------------------- step 4 --------------------------------- #

  B <- list()
  B[[1]] <- B0hat
  for (i in 1:nlags) {
    B[[i+1]] <- (diag(dim(B0hat)[1]) - B0hat) %*% VARres$Mhat[[i]]
  }


  # -------------------------------- return --------------------------------- #

  res <- list()
  res$resid <- VARres$residuals
  res$const <- VARres$const
  res$Mhat <- VARres$Mhat
  res$Bhat <- B
  res$var_order <- reslg$k

  res

}
