main_FirmGrowth <- function(path, nlags=1) {

  # This code reproduces the results of Chapter 4 of the OBES paper
  #    'Causal Inference by Independent Component Analysis:
  #                Theory and Applications'
  # by A. Moneta, D. Entner, P.O. Hoyer, and A. Coad

  # nlags: integer, number of time lags included in the VAR/SVAR model
  #        in paper either 1 or 2
  # path: path to the data file including file name, f.ex. "myPath/theFile.txt"

  # ------ Load Data ------
  # the data can be bought from the Compustat database
  # preprocessing as stated in the paper
  # columns: Firmname, Year, Empl.gr, Sales.gr, R&Dexp.gr, Profit.gr
  Growth <- read.table(path, colClasses=c('character', 'character', 'numeric',
            'real', 'real', 'real', 'real'))

  # ------ Correlation between the variables ------
  corr <- array(0,dim=c(4,4))
  pval <- array(0,dim=c(4,4))
  for (i in 3:6) {
    for (j in i:6) {
      temp <- unclass(cor.test(Growth[,i],Growth[,j],method=c("pearson")))
      corr[i-2,j-2] <- temp$estimate
      pval[i-2,j-2] <- temp$p.value
    }
  }

  print("Correlation between growth variables:")
  print(corr)

  # ------ Estimate VAR-LiNGAM ------
  Growth_can <- canonicalform_USfirmsdata(Growth,nlags=nlags)

  result <- VARLiNGAM(Growth_can,"lad",ntests=TRUE)

  cat('\n ------------------------------------------------------- \n')
  cat('VAR estimates (A1, ..., Ap): \n')
  print(result$Mhat)

  cat('VAR-LiNGAM estimates (Gamma0, Gamma1, ..., Gammap): \n')
  print(result$Bhat)

  nodelabels <- list()
  nodelabels[[1]] <- "Empl.gr"
  nodelabels[[2]] <- "Sales.gr"
  nodelabels[[3]] <- "RnD.gr"
  nodelabels[[4]] <- "Opinc.gr"

  # ------ Bootstrap ------
  bootstrapping = TRUE
  if (bootstrapping) {
    # notation in paper: B = B_0, Gamma_i = B_i, A_i=M_i, i>0
    cat('\n ------------------------------------------------------- \n')
    cat('Start bootstrap analysis. \n')

    nboot <- 100

    if (nlags==1) {
      mySeed <- 4008
    }
    else {
      if (nlags==2) {
        mySeed <- 6293
      }
      else{
        mySeed <- ceiling(runif(1)*10000) 
      }
    }
    
    set.seed(mySeed)
    boot_sd_firmgrowth(Growth_can, nboot, result$Mhat, result$Bhat)
    cat("used seed :", mySeed, "\n")

  } # end bootstrapping

}
