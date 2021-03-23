# lingam - perform the LiNGAM analysis
#
# SYNTAX:
# res <- lingam( X )
#
# INPUT:
# X     - Data matrix: each row is an observed variable, each
#         column one observed sample. The number of columns should
#         be far larger than the number of rows.
#
# OUTPUT:
# res$B     - Matrix of estimated connection strenghts
# res$stde  - Standard deviations of disturbance variables
# res$ci    - Constants
# res$k     - An estimated causal ordering
#
# (Note that B, stde, and ci are all ordered according to the same
#  variable ordering as X.)
#
# Version: 0.1 (8 Jan 2008) - first version for R
#
# See also ESTIMATE, PRUNE.

lingam <- function( X ) {

  temp <- estimate( X )
#   print(temp$k)
#   k <- c(2,1,4,3,6,5) for macro-data application
#   print(k)
  res <- prune( X, temp$k )

  res$k <- temp$k
#   res$W <- temp$W
#   res$ICs <- temp$ICs

  res

#   res1 <- list()
#   res1 <- list(Bpruned=res$Bpruned, stde=res$stde, ci=res$ci, W = temp$W,
#                ICs = temp$ICs, k=temp$k)
#   res1
  
}
