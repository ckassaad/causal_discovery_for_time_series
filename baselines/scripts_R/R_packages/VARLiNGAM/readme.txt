Codepackage to the paper by A. Moneta, D. Entner, P.O. Hoyer, and A. Coad

         'Causal Inference by Independent Component Analysis:
                     Theory and Applications'
      


The code is written in R, see http://cran.r-project.org/


Required packages:

nortest   (for Gauss_Tests.R)
tseries   (for Gauss_Tests.R)
fastICA   (for estimate.R (LiNGAM))
quantreg  (for lad_est.R)
sm        (for fm_LAD.R)

To load the required files, switch to the codepackage directory and type

  source("sourcedir.R")
  source("main1.R")


1. Replicating the results of the paper

   By calling the main1 function, the results can be dublicated as follows. Make
   sure to have the right path to the datafile (in txt format) in the variable
   path (f.ex. path <- "myPath/datafile.txt")

   Section 4 - Firmgwoth:
     main1('micro', path, nlags=1) # for the 1-lag analysis
     main1('micro', path, nlags=2) # for the 2-laga analysis

   Section 5 - Monetary Policy:
     # this includes all the analysis of the paper
     main1('macro', path, boot_sd=TRUE, boot_irf=TRUE, subsamples=TRUE)
     # alternatively, to save some time and only get results without bootstrap
     main1('macro', path, boot_sd=FALSE, boot_irf=FALSE, subsamples=TRUE)


2. Applying VAR-LiNGAM to own data

   Load your data into R under the name X (f.ex. X <- read.table('myDataFile'))
   with the variables in the columns and the observations in the rows. 
   Call the functions (with nlgas being the number lags to be included in the
   analysis)

     main1('own') # load the files
     X_can <- tsdata2canonicalform(X,nlags) # put data into canonical form
     result <- VARLiNGAM(X_can,"ols", pruning=FALSE)

   In the variable result are the VAR coefficients, the VAR-LiNGAM coefficients
   the residuals, and the estiamted causal order. See the file VARLiNGAM.R for
   details on input and output.


For questions please contact Doris Entner (doris.entner'at'cs.helsinki.fi).

