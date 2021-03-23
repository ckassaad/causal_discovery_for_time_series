lad_est <- function(X,y,meth="qreg",regstats=FALSE) {

  # least absolute deviation estimation of a regression model
  #
  # INPUT
  # X, y: Data such that of y=c+Xb (with c, b the parameters to be estimated)
  # meth: "qreg" using quantile regression (needs package 'quantreg')
  #       "iter" using iterative algorithm for LAD estimation, see for example 
  #              2004 MPRA: LAD estimation, p.6
  # regstats: boolean, if TRUE print t-statistics, standard errors, p-values
  # 
  # OUTPUT
  # [c b]: regression coefficients s.t. y=c+Xb

  if (meth == "qreg") {

    require(quantreg)

    dy <- dim(y)
    bhat <- array(0, dim=c(dy[2],dim(X)[2]+1) )
  
    for (i in 1:dy[2]) { # for each equation separately
      quanreg <- rq(y[,i]~X, tau=0.5, method="br")#"fn")
      bhat[i,] <- quanreg$coeff

      if (regstats) {
        print(summary(quanreg,se="iid")) # Reg-Coeffs, Std. Errors and t-values
                                # iid, nid, ker, boot
        yhat <- bhat[i,1]+X%*%bhat[i,2:(dim(X)[2]+1)]
     
        # R2 - goodness-of-fit
        R2 <- 1 - sum(abs(y[,i]-yhat)) / sum(abs(y[,i]-median(y[,i])))
        cat("R2 =", R2, "\n")
      }

    }

    return(bhat)
  }


  if (meth == "iter") {

    X <- cbind(matrix(1,nrow=dim(X)[1],ncol=1),X)

    bhat <- solve(t(X)%*%X)%*%t(X)%*%y # OLS-estimator
    yhat <- X%*%bhat
    ehat <- y-yhat
    bhat_old <- bhat+1 # that it jumps in while-loop
    eps <- 0.005 #0.01
    dims <- dim(X)

    max_loop <- 1000

    for (i in 1:dim(y)[2]) { # for each equation separately

      cnt <- 0

      while ( sum(abs(bhat[,i]-bhat_old[,i])) > eps && cnt < max_loop) {

        bool <- abs(ehat[,i]) < 0.00001
        tmp <- ehat[,i]
        # if residuals are too small set to 0.00001 (use approach of Fair)
        tmp[bool] <- 0.00001
        # elementwise division of each column = t(X)*W
        Temp <- t(X/as.vector(abs(tmp)))

        bhat_old[,i] <- bhat[,i]
        bhat[,i] <- solve(Temp%*%X)%*%Temp%*%y[,i]
        yhat[,i] <- X%*%bhat[,i]
        ehat[,i] <- y[,i]-yhat[,i]
        cnt <- cnt+1

        if (cnt == max_loop) {
          cat('WARNING: not yet under threshold of', eps, 'for last two consecutive iterations. Value is', sum(abs(bhat[,i]-bhat_old[,i])), 'but maximal number of iterations' , max_loop,'reached.\n')
          cat('Problem occured for variable ', i, '.\n')
        }
      }
    }

    return(t(bhat))
  }

}
