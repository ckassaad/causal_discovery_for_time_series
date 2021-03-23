# lingamer:
# Runs LiNGAM on the data and returns the output as a list of DAGs
# and their corresponding scores, in the same format as bayeslingam
# and dealer, for easy comparison
#
# SYNTAX:
# r <- lingamer( X, dags=NULL )
#
# INPUT:
# X         - data matrix (each row is one data vector, 20 samples minimum!)
# dags      - list of dags (default: output full list of all DAGs)
#
# OUTPUT:
# r$DAGs    - contains the DAGs
# r$loglike - log-likelihood values of DAGs (0 or infinite)
# r$prob    - probabilities of DAGs (1 or 0)
#

lingamer <- function( X ) {

  # number of variables
  nvars<-ncol(X)


  # initialize
  r <- list()
  r$loglike<-0
  r$prob <- 1

  # if there are not enough samples for lingam to run without crashing
  # then repeat the data (a bit of a hack, but works)  
  Xoriginal<-X
  samples_needed<-10+nvars*10

  while ( nrow(X) < samples_needed ) {
    X<-rbind(X,Xoriginal)
  }
 # print(X)
  # run lingam
  R<-lingam(t(X))

  # threshold (we care only about the structure here)
  B<- (abs(R$Bprune) > 1e-5)*1

  r$DAGs<-bdag.to.cdag(B)

  # set (loglike, prob) of correct DAG to (0,1), others (-Inf,0)
#  for (i in 1:nrow(dags) ) {
#    Bdag<-dagvector2B(dags[i,])
#    if ( all(Bdag == B ) ) {
#      r$loglike[i]<-0
#      r$prob[i]<-1
#    } else {
#      r$loglike[i]<-(-Inf)
#      r$prob[i]<-0
#    }
#  }

  # return the result
  r
}
