macro_subsamples <- function(Data, est_meth, fm_LAD=FALSE) {

  # Robustness analysis, Section 5.3

  # INPUT
  # Data: origianl Data set with variables GDP, PGDP, NBR, BR, FFR, PSCCOM
  # est_meth: one of "ols", "lad", "vecm"
  # fm_LAD: boolean, when est_meth="lad" and fm_LAD=TRUE use fully modified LAD
  #
  # OUTPUT
  # prints matrices of instantaneous effects of different subsamples:
  #     whole sample
  #     1965:1 - 1979:9
  #     1979:10 -1996:12
  #     1984:2 - 1996:12
  #     1988:9 - 1996:12
  # plots impulse response functions as in Figure 8

  # subsamples
  Data1 <- Data[1:177,]            # 1965:1 - 1979:9
  Data2 <- Data[178:nrow(Data),]   # 1979:10 -1996:12
  Data3 <- Data[230:nrow(Data),]   # 1984:2 - 1996:12
  Data4 <- Data[285:nrow(Data),]   # 1988:9 - 1996:12

  # put data into canonical form with 7 lags
  nlags <- 7
  Data_can <- tsdata2canonicalform(Data,nlags)
  Data_can1 <- tsdata2canonicalform(Data1,nlags)
  Data_can2 <- tsdata2canonicalform(Data2,nlags)
  Data_can3 <- tsdata2canonicalform(Data3,nlags)
  Data_can4 <- tsdata2canonicalform(Data4,nlags)

  # Apply VAR-LiNGAM to the data sets
  result <- VARLiNGAM(Data_can, est_meth, ntests=FALSE, corank=3, fmlad=fm_LAD,
            pruning=FALSE)
  result1 <- VARLiNGAM(Data_can1, est_meth, ntests=FALSE, corank=3,
             fmlad=fm_LAD, pruning=FALSE)
  result2 <- VARLiNGAM(Data_can2, est_meth, ntests=FALSE, corank=3,
             fmlad=fm_LAD, pruning=FALSE)
  result3 <- VARLiNGAM(Data_can3, est_meth, ntests=FALSE, corank=3,
             fmlad=fm_LAD, pruning=FALSE)
  result4 <- VARLiNGAM(Data_can4, est_meth, ntests=FALSE, corank=3,
             fmlad=fm_LAD, pruning=FALSE)

  # reduced form residuals
  u_res <- result$resid 
  u_res1 <- result1$resid 
  u_res2 <- result2$resid 
  u_res3 <- result3$resid 
  u_res4 <- result4$resid 

  nodelabels <- colnames(Data)

  # instantaneous effect matrices
  B0hat <- result$Bhat[[1]]
  ord <- c(2,1,4,3,6,5)
  B0 <- choleski(u_res,ord)
  print("Results for whole sample: causal order, B0hat, B0hat with original order")
  print(nodelabels[result$var_order])
  print(round(B0hat,4))
  print(round(B0,4))

  B0hat1 <- result1$Bhat[[1]]
  ord1 <- c(2,1,4,3,6,5)
  B01 <- choleski(u_res1,ord)
  print("Results for 1965:1 - 1979:9: causal order, B0hat, B0hat with original order")
  print(nodelabels[result1$var_order])
  print(round(B0hat1,4))
  print(round(B01,4))

  B0hat2 <- result2$Bhat[[1]]
  ord2 <- c(2,1,4,3,6,5)
  B02 <- choleski(u_res2,ord)
  print("Results for 1979:10 -1996:12: causal order, B0hat, B0hat with original order")
  print(nodelabels[result2$var_order])
  print(round(B0hat2,4))
  print(round(B02,4))

  B0hat3 <- result3$Bhat[[1]]
  ord3 <- c(2,1,4,3,6,5)
  B03 <- choleski(u_res3,ord)
  print("Results for 1984:2 - 1996:12: causal order, B0hat, B0hat with original order")
  print(nodelabels[result3$var_order])
  print(round(B0hat3,4))
  print(round(B03,4))

  B0hat4 <- result4$Bhat[[1]]
  ord4 <- c(2,1,4,3,6,5)
  B04 <- choleski(u_res4,ord)
  print("Results for 1988:9 - 1996:12: causal order, B0hat, B0hat with original order")
  print(nodelabels[result4$var_order])
  print(round(B0hat4,4))
  print(round(B04,4))

  # coefficient matrices of reduced form VAR
  AA <- result$Mhat
  AA1 <- result1$Mhat
  AA2 <- result2$Mhat
  AA3 <- result3$Mhat
  AA4 <- result4$Mhat


  # --------------------------------------------------------------------------
  # Impulse Response Functions

  lag_irf <- 25
  IRF <- irf(AA, B0, lag_irf, u_res)
  IRF1 <- irf(AA1, B01, lag_irf, u_res1)
  IRF2 <- irf(AA2, B02, lag_irf, u_res2)
  IRF3 <- irf(AA3, B03, lag_irf, u_res3)
  IRF4 <- irf(AA4, B04, lag_irf, u_res4)

  x11(width=12, height=12)
  par(mfrow=c(3,3))
  nd <- colnames(Data)
  n <- 0
  for (i in 1:6){
    for (j in 1:6){

      n <- n+1   

      if (is.element(i,c(1,2,5)) && is.element(j,c(3,4,5)) ) {

	cirf <- c(IRF[n,],IRF1[n,], IRF2[n,], IRF3[n,], IRF4[n,])
	di <- abs(max(cirf) - min(cirf))
	miny <- min(cirf) - (di/2)
	maxy <- max(cirf) + (di/2)

	matplot(0:25, IRF[n,], t="l", ylim=c(miny,maxy), xlab="", ylab="",
                lty=1)
	matplot(0:25, IRF1[n,], t="l", col=2, lty=2, add=TRUE)
	matplot(0:25, IRF2[n,], t="l", col=3, lty=3, add=TRUE)
	matplot(0:25, IRF3[n,], t="l", col=4, lty=4, add=TRUE)
	matplot(0:25, IRF4[n,], t="l", col=5, lty=5, add=TRUE)
	abline(h=0, col="grey", lwd=0.5)
	legend(12, maxy+ maxy/15, c("S0: 1965:1 - 1996:12", 
	      "S1: 1965:1 - 1979:9","S2: 1979:10-1996:12",
	      "S3: 1984:2 - 1996:12", "S4: 1988:9 - 1996:12"),
	      lty=c(1,2,3,4,5), col=c(1,2,3,4,5), bty="n") 
	title(paste("Responses of", nd[i], "to", nd[j], "(subsamples)") )

      } # end if

    }
  }

}
