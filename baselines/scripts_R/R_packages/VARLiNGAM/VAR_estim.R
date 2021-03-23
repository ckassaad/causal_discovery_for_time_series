VAR_estim <- function(Data, est_meth="ols", regstats=TRUE, corank=NA,
             fmlad=FALSE) {

  # estimation of a VAR model
  #
  # INPUT
  # Data: time series data in canonical form
  #       (see function 'tsdata2canonicalform.R')
  # est_meth: string giving the estimation method of the levels:
  #       "ols" for for ordinary least squares
  #       "lad" for least absolute deviation
  #       "vecm" for an Vector Error Correction Model
  # regstats: return regression statistics, for est_meth="ols" or "lad" only
  # corank: int, only considered when est_meth="vecm", give the cointegraion
  #        rank in the VECM
  # fmlad: boolen, only considered when est_meth="lad"
  #        if TRUE the fully modified LAD is used (Phillips ET 1995)
  #
  # OUTPUT
  # List with estimated constant term and coefficients for regression and the
  #      residuals:
  #   $const (constant term in VAR)
  #   $Mhat (list of length nlags of estimated coeff. matrices for VAR)
  #   $residuals (estimated residuals)

  dims <- dim(Data)
  nobs <- dims[1] # number of observations
  nvar <- dims[2] - which(colnames(Data)=="curval1") + 1 # number of variables
  nlags <- (dims[2] - 2)/nvar - 1 # time lags included

  # get data such that y=c+Xb (with c, b the parameters to be estimated)
  X <- Data[,3:(dims[2]-nvar)]
  y <- Data[,(dims[2]-nvar+1):dims[2]]

  # ------ Estimation ------

  if (est_meth == "ols") {
    est <- ols_est(X,y,regstats) # ols-estimates:, nvar x (1+nvar*nlags) matrix

    # save results
    constant <- as.matrix(est[,1])
    Mhat <- list()
    yhat <- array(constant,dim=c(nvar,nobs))
    for (i in 1:nlags) {
      Mhat[[i]] <- as.matrix(est[,(2+(i-1)*nvar):(1+i*nvar)])
      yhat <- yhat + Mhat[[i]]%*%t(X[,(1+(i-1)*nvar):(i*nvar)])
    }
  }

  if (est_meth == "lad") {
    est <- lad_est(X,y,"qreg",regstats) # lad estimates

    # save results
    constant <- as.matrix(est[,1])
    Mhat <- list()
    yhat <- array(constant,dim=c(nvar,nobs))
    for (i in 1:nlags) {
      Mhat[[i]] <- as.matrix(est[,(2+(i-1)*nvar):(1+i*nvar)])
      yhat <- yhat + Mhat[[i]]%*%t(X[,(1+(i-1)*nvar):(i*nvar)])
    }

    if (fmlad) { # if we use fully modified LAD
      print("Getting fully modified LAD estimates. This takes some time.")
      res <- y - t(yhat) # residuals from LAD estimation
      yhat <- array(constant,dim=c(nvar,nobs))
      # Data in canonical form... put it back in 'normal' Form
      Dn <- rbind(X[1:nlags,(dim(X)[2]-nvar+1):dim(X)[2]],y)
      fmlad <- fm_LAD(Dn, res) # calculate residuals with fm-lad
      for (i in 1:nlags) {
        Mhat[[i]] <- Mhat[[i]] + fmlad[[i]]
        yhat <- yhat + Mhat[[i]]%*%t(X[,(1+(i-1)*nvar):(i*nvar)])
      }
      print("Done.")
    }

  }

  if (est_meth=="vecm") { # VECM estimation

      # Data in canonical form... put it back in 'normal' Form
      Dn <- rbind(X[1:nlags,(dim(X)[2]-nvar+1):dim(X)[2]],y)
      Mhat <- list()

      vecm <- co_ml(Dn, nlags, corank)
      constant <- vecm$const
      resids <- vecm$residuals
      Mhat<-vecm$AA
      yhat <- array(constant,dim=c(nvar,nobs)) # calculate residuals 
      for (i in 1:nlags) {
        yhat <- yhat + Mhat[[i]]%*%t(X[,(1+(i-1)*nvar):(i*nvar)])
      }
  }


  ls <- list()
  ls$const <- constant
  ls$Mhat <- Mhat # estimated coefficient matrices for VAR
  ls$residuals <- y - t(yhat)
  ls

}
