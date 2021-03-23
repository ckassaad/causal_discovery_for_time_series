# This file contains the commands to run the simulations on discrete data
# in "On Causal discovery from Time-Series Data using tsFCI"
# 
# Copyright (C) 2010 Doris Entner
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see http://www.gnu.org/licenses/



# --------------------------------------------------------------------------- #
# ----------------------- simulations on discrete data ---------------------- #
# --------------------------------------------------------------------------- #

simDataDisc <- function(q=0.25, alg="tscfci",nrep=4, sig=0.01, nsim=100) {

  # nrep = tau+1 in paper

  path <- "./simulation_results/" # where results will be saved
  sig_write <- sig*100

  if (q!=0.25 & q!=0.5) {
    return("Input parameter q must either be 0.25 or 0.5.")
  }

# ------------------------- edgepercentage q = 0.25 ------------------------- #
  if (q==0.25) {

    # this does "nsim" simulations for each setting
    # saves generating models and resulting pags for each setting in a file in 
    #   the folder simulation_results
    # also saves a summary of the socres in the folder simulation_results
    # can load files using readListFromFile(filenname) - see main_tetrad_fci.R

    seeds <- c(7834,2175,6273,5940,936,8722,9355,195)
    n <- length(seeds) # number of different seetings for nobs and nhid

    # --- 100 samples

    simOut21 <- doSimulations_universal(2,1,0.25,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[1])
    simOut31 <- doSimulations_universal(3,1,0.25,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[2])
    simOut32 <- doSimulations_universal(3,2,0.25,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[3])
    simOut42 <- doSimulations_universal(4,2,0.25,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[4])
    simOut43 <- doSimulations_universal(4,3,0.25,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[5])
    simOut52 <- doSimulations_universal(5,2,0.25,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[6])
    simOut53 <- doSimulations_universal(5,3,0.25,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[7])
    simOut63 <- doSimulations_universal(6,3,0.25,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[8])

    scores21 <- checkSimulationOutput_scores(simOut21,2,nrep)$scores
    scores31 <- checkSimulationOutput_scores(simOut31,3,nrep)$scores
    scores32 <- checkSimulationOutput_scores(simOut32,3,nrep)$scores
    scores42 <- checkSimulationOutput_scores(simOut42,4,nrep)$scores
    scores43 <- checkSimulationOutput_scores(simOut43,4,nrep)$scores
    scores52 <- checkSimulationOutput_scores(simOut52,5,nrep)$scores
    scores53 <- checkSimulationOutput_scores(simOut53,5,nrep)$scores
    scores63 <- checkSimulationOutput_scores(simOut63,6,nrep)$scores

    # direct-cause scores
    scores100DC <- array(0,dim=c(n,2))
    scores100DC[1,] <- scores21[1,]
    scores100DC[2,] <- scores31[1,]
    scores100DC[3,] <- scores32[1,]
    scores100DC[4,] <- scores42[1,]
    scores100DC[5,] <- scores43[1,]
    scores100DC[6,] <- scores52[1,]
    scores100DC[7,] <- scores53[1,]
    scores100DC[8,] <- scores63[1,]

    # ancestor scores
    scores100AN <- array(0,dim=c(n,2))
    scores100AN[1,] <- scores21[2,]
    scores100AN[2,] <- scores31[2,]
    scores100AN[3,] <- scores32[2,]
    scores100AN[4,] <- scores42[2,]
    scores100AN[5,] <- scores43[2,]
    scores100AN[6,] <- scores52[2,]
    scores100AN[7,] <- scores53[2,]
    scores100AN[8,] <- scores63[2,]

    # pairwise scores
    scores100PW <- array(0,dim=c(n,2))
    scores100PW[1,] <- scores21[3,]
    scores100PW[2,] <- scores31[3,]
    scores100PW[3,] <- scores32[3,]
    scores100PW[4,] <- scores42[3,]
    scores100PW[5,] <- scores43[3,]
    scores100PW[6,] <- scores52[3,]
    scores100PW[7,] <- scores53[3,]
    scores100PW[8,] <- scores63[3,]

    # marginalize over number of observed and hidden variables
    score100 <- array(0,dim=c(3,2))
    score100[1,] <- colSums(scores100DC)/n
    score100[2,] <- colSums(scores100AN)/n
    score100[3,] <- colSums(scores100PW)/n

    # --- 1000 samples

    simOut21 <- doSimulations_universal(2,1,0.25,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[1])
    simOut31 <- doSimulations_universal(3,1,0.25,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[2])
    simOut32 <- doSimulations_universal(3,2,0.25,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[3])
    simOut42 <- doSimulations_universal(4,2,0.25,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[4])
    simOut43 <- doSimulations_universal(4,3,0.25,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[5])
    simOut52 <- doSimulations_universal(5,2,0.25,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[6])
    simOut53 <- doSimulations_universal(5,3,0.25,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[7])
    simOut63 <- doSimulations_universal(6,3,0.25,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[8])

    scores21 <- checkSimulationOutput_scores(simOut21,2,nrep)$scores
    scores31 <- checkSimulationOutput_scores(simOut31,3,nrep)$scores
    scores32 <- checkSimulationOutput_scores(simOut32,3,nrep)$scores
    scores42 <- checkSimulationOutput_scores(simOut42,4,nrep)$scores
    scores43 <- checkSimulationOutput_scores(simOut43,4,nrep)$scores
    scores52 <- checkSimulationOutput_scores(simOut52,5,nrep)$scores
    scores53 <- checkSimulationOutput_scores(simOut53,5,nrep)$scores
    scores63 <- checkSimulationOutput_scores(simOut63,6,nrep)$scores

    # direct-cause scores
    scores1000DC <- array(0,dim=c(n,2))
    scores1000DC[1,] <- scores21[1,]
    scores1000DC[2,] <- scores31[1,]
    scores1000DC[3,] <- scores32[1,]
    scores1000DC[4,] <- scores42[1,]
    scores1000DC[5,] <- scores43[1,]
    scores1000DC[6,] <- scores52[1,]
    scores1000DC[7,] <- scores53[1,]
    scores1000DC[8,] <- scores63[1,]

    # ancestor scores
    scores1000AN <- array(0,dim=c(n,2))
    scores1000AN[1,] <- scores21[2,]
    scores1000AN[2,] <- scores31[2,]
    scores1000AN[3,] <- scores32[2,]
    scores1000AN[4,] <- scores42[2,]
    scores1000AN[5,] <- scores43[2,]
    scores1000AN[6,] <- scores52[2,]
    scores1000AN[7,] <- scores53[2,]
    scores1000AN[8,] <- scores63[2,]

    # pairwise scores
    scores1000PW <- array(0,dim=c(n,2))
    scores1000PW[1,] <- scores21[3,]
    scores1000PW[2,] <- scores31[3,]
    scores1000PW[3,] <- scores32[3,]
    scores1000PW[4,] <- scores42[3,]
    scores1000PW[5,] <- scores43[3,]
    scores1000PW[6,] <- scores52[3,]
    scores1000PW[7,] <- scores53[3,]
    scores1000PW[8,] <- scores63[3,]

    # marginalize over number of observed and hidden variables
    score1000 <- array(0,dim=c(3,2))
    score1000[1,] <- colSums(scores1000DC)/n
    score1000[2,] <- colSums(scores1000AN)/n
    score1000[3,] <- colSums(scores1000PW)/n

    # --- 10000 samples

    simOut21 <- doSimulations_universal(2,1,0.25,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[1])
    simOut31 <- doSimulations_universal(3,1,0.25,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[2])
    simOUt32 <- doSimulations_universal(3,2,0.25,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[3])
    simOut42 <- doSimulations_universal(4,2,0.25,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[4])
    simOut43 <- doSimulations_universal(4,3,0.25,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[5])
    simOut52 <- doSimulations_universal(5,2,0.25,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[6])
    simOut53 <- doSimulations_universal(5,3,0.25,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[7])
    simOut63 <- doSimulations_universal(6,3,0.25,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[8])

    scores21 <- checkSimulationOutput_scores(simOut21,2,nrep)$scores
    scores31 <- checkSimulationOutput_scores(simOut31,3,nrep)$scores
    scores32 <- checkSimulationOutput_scores(simOut32,3,nrep)$scores
    scores42 <- checkSimulationOutput_scores(simOut42,4,nrep)$scores
    scores43 <- checkSimulationOutput_scores(simOut43,4,nrep)$scores
    scores52 <- checkSimulationOutput_scores(simOut52,5,nrep)$scores
    scores53 <- checkSimulationOutput_scores(simOut53,5,nrep)$scores
    scores63 <- checkSimulationOutput_scores(simOut63,6,nrep)$scores

    # direct-cause scores
    scores10000DC <- array(0,dim=c(n,2))
    scores10000DC[1,] <- scores21[1,]
    scores10000DC[2,] <- scores31[1,]
    scores10000DC[3,] <- scores32[1,]
    scores10000DC[4,] <- scores42[1,]
    scores10000DC[5,] <- scores43[1,]
    scores10000DC[6,] <- scores52[1,]
    scores10000DC[7,] <- scores53[1,]
    scores10000DC[8,] <- scores63[1,]

    # ancestor scores
    scores10000AN <- array(0,dim=c(n,2))
    scores10000AN[1,] <- scores21[2,]
    scores10000AN[2,] <- scores31[2,]
    scores10000AN[3,] <- scores32[2,]
    scores10000AN[4,] <- scores42[2,]
    scores10000AN[5,] <- scores43[2,]
    scores10000AN[6,] <- scores52[2,]
    scores10000AN[7,] <- scores53[2,]
    scores10000AN[8,] <- scores63[2,]

    # pairwise scores
    scores10000PW <- array(0,dim=c(n,2))
    scores10000PW[1,] <- scores21[3,]
    scores10000PW[2,] <- scores31[3,]
    scores10000PW[3,] <- scores32[3,]
    scores10000PW[4,] <- scores42[3,]
    scores10000PW[5,] <- scores43[3,]
    scores10000PW[6,] <- scores52[3,]
    scores10000PW[7,] <- scores53[3,]
    scores10000PW[8,] <- scores63[3,]

    # marginalize over number of observed and hidden variables
    score10000 <- array(0,dim=c(3,2))
    score10000[1,] <- colSums(scores10000DC)/n
    score10000[2,] <- colSums(scores10000AN)/n
    score10000[3,] <- colSums(scores10000PW)/n

    # --- write to file, before marginalization
    res <- list(s100DC=scores100DC,s1000DC=scores1000DC,s10000DC=scores10000DC,
                s100AN=scores100AN,s1000AN=scores1000AN,s10000AN=scores10000AN,
                s100PW=scores100PW,s1000PW=scores1000PW,s10000PW=scores10000PW)
    q_write <- q*100
    filename <- paste(path,"RESdisc",q_write,alg,nrep,"_",sig_write,".txt",
      sep="")
    writeListToFile(filename,res)

    return(list(s100=score100, s1000=score1000, s10000=score10000))

  } # end q=0.25

# ------------------------- edgepercentage q = 0.50 ------------------------- #
  if (q==0.5) {

    # this does "nsim" simulations for each setting
    # saves generating models and resulting pags for each setting in a file in 
    #   the folder simulation_results
    # also saves a summary of the socres in the folder simulation_results
    # can load files using readListFromFile(filenname) - see main_tetrad_fci.R

    seeds <- c(2113,7526,6642,7173,3035,5623)
    n <- length(seeds) # number of different seetings for nobs and nhid

    # --- 100 samples

    simOut21 <- doSimulations_universal(2,1,0.5,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[1])
    simOut31 <- doSimulations_universal(3,1,0.5,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[2])
    simOut32 <- doSimulations_universal(3,2,0.5,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[3])
    simOut42 <- doSimulations_universal(4,2,0.5,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[4])
    simOut43 <- doSimulations_universal(4,3,0.5,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[5])
    simOut52 <- doSimulations_universal(5,2,0.5,alg,nrep,nsim,"discrete",100,
                sig,FALSE,seeds[6])

    scores21 <- checkSimulationOutput_scores(simOut21,2,nrep)$scores
    scores31 <- checkSimulationOutput_scores(simOut31,3,nrep)$scores
    scores32 <- checkSimulationOutput_scores(simOut32,3,nrep)$scores
    scores42 <- checkSimulationOutput_scores(simOut42,4,nrep)$scores
    scores43 <- checkSimulationOutput_scores(simOut43,4,nrep)$scores
    scores52 <- checkSimulationOutput_scores(simOut52,5,nrep)$scores

    # direct-cause scores
    scores100DC <- array(0,dim=c(n,2))
    scores100DC[1,] <- scores21[1,]
    scores100DC[2,] <- scores31[1,]
    scores100DC[3,] <- scores32[1,]
    scores100DC[4,] <- scores42[1,]
    scores100DC[5,] <- scores43[1,]
    scores100DC[6,] <- scores52[1,]

    # ancestor scores
    scores100AN <- array(0,dim=c(n,2))
    scores100AN[1,] <- scores21[2,]
    scores100AN[2,] <- scores31[2,]
    scores100AN[3,] <- scores32[2,]
    scores100AN[4,] <- scores42[2,]
    scores100AN[5,] <- scores43[2,]
    scores100AN[6,] <- scores52[2,]

    # pairwise scores
    scores100PW <- array(0,dim=c(n,2))
    scores100PW[1,] <- scores21[3,]
    scores100PW[2,] <- scores31[3,]
    scores100PW[3,] <- scores32[3,]
    scores100PW[4,] <- scores42[3,]
    scores100PW[5,] <- scores43[3,]
    scores100PW[6,] <- scores52[3,]

    # marginalize over number of observed and hidden variables
    score100 <- array(0,dim=c(3,2))
    score100[1,] <- colSums(scores100DC)/n
    score100[2,] <- colSums(scores100AN)/n
    score100[3,] <- colSums(scores100PW)/n

    # --- 1000 samples

    simOut21 <- doSimulations_universal(2,1,0.5,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[1])
    simOut31 <- doSimulations_universal(3,1,0.5,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[2])
    simOut32 <- doSimulations_universal(3,2,0.5,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[3])
    simOut42 <- doSimulations_universal(4,2,0.5,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[4])
    simOut43 <- doSimulations_universal(4,3,0.5,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[5])
    simOut52 <- doSimulations_universal(5,2,0.5,alg,nrep,nsim,"discrete",1000,
                sig,FALSE,seeds[6])

    scores21 <- checkSimulationOutput_scores(simOut21,2,nrep)$scores
    scores31 <- checkSimulationOutput_scores(simOut31,3,nrep)$scores
    scores32 <- checkSimulationOutput_scores(simOut32,3,nrep)$scores
    scores42 <- checkSimulationOutput_scores(simOut42,4,nrep)$scores
    scores43 <- checkSimulationOutput_scores(simOut43,4,nrep)$scores
    scores52 <- checkSimulationOutput_scores(simOut52,5,nrep)$scores

    # direct-cause scores
    scores1000DC <- array(0,dim=c(n,2))
    scores1000DC[1,] <- scores21[1,]
    scores1000DC[2,] <- scores31[1,]
    scores1000DC[3,] <- scores32[1,]
    scores1000DC[4,] <- scores42[1,]
    scores1000DC[5,] <- scores43[1,]
    scores1000DC[6,] <- scores52[1,]

    # ancestor scores
    scores1000AN <- array(0,dim=c(n,2))
    scores1000AN[1,] <- scores21[2,]
    scores1000AN[2,] <- scores31[2,]
    scores1000AN[3,] <- scores32[2,]
    scores1000AN[4,] <- scores42[2,]
    scores1000AN[5,] <- scores43[2,]
    scores1000AN[6,] <- scores52[2,]

    # pairwise scores
    scores1000PW <- array(0,dim=c(n,2))
    scores1000PW[1,] <- scores21[3,]
    scores1000PW[2,] <- scores31[3,]
    scores1000PW[3,] <- scores32[3,]
    scores1000PW[4,] <- scores42[3,]
    scores1000PW[5,] <- scores43[3,]
    scores1000PW[6,] <- scores52[3,]

    # marginalize over number of observed and hidden variables
    score1000 <- array(0,dim=c(3,2))
    score1000[1,] <- colSums(scores1000DC)/n
    score1000[2,] <- colSums(scores1000AN)/n
    score1000[3,] <- colSums(scores1000PW)/n

    # --- 10000 samples

    simOut21 <- doSimulations_universal(2,1,0.5,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[1])
    simOut31 <- doSimulations_universal(3,1,0.5,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[2])
    simOUt32 <- doSimulations_universal(3,2,0.5,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[3])
    simOut42 <- doSimulations_universal(4,2,0.5,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[4])
    simOut43 <- doSimulations_universal(4,3,0.5,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[5])
    simOut52 <- doSimulations_universal(5,2,0.5,alg,nrep,nsim,"discrete",10000,
                sig,FALSE,seeds[6])

    scores21 <- checkSimulationOutput_scores(simOut21,2,nrep)$scores
    scores31 <- checkSimulationOutput_scores(simOut31,3,nrep)$scores
    scores32 <- checkSimulationOutput_scores(simOut32,3,nrep)$scores
    scores42 <- checkSimulationOutput_scores(simOut42,4,nrep)$scores
    scores43 <- checkSimulationOutput_scores(simOut43,4,nrep)$scores
    scores52 <- checkSimulationOutput_scores(simOut52,5,nrep)$scores

    # direct-cause scores
    scores10000DC <- array(0,dim=c(n,2))
    scores10000DC[1,] <- scores21[1,]
    scores10000DC[2,] <- scores31[1,]
    scores10000DC[3,] <- scores32[1,]
    scores10000DC[4,] <- scores42[1,]
    scores10000DC[5,] <- scores43[1,]
    scores10000DC[6,] <- scores52[1,]

    # ancestor scores
    scores10000AN <- array(0,dim=c(n,2))
    scores10000AN[1,] <- scores21[2,]
    scores10000AN[2,] <- scores31[2,]
    scores10000AN[3,] <- scores32[2,]
    scores10000AN[4,] <- scores42[2,]
    scores10000AN[5,] <- scores43[2,]
    scores10000AN[6,] <- scores52[2,]

    # pairwise scores
    scores10000PW <- array(0,dim=c(n,2))
    scores10000PW[1,] <- scores21[3,]
    scores10000PW[2,] <- scores31[3,]
    scores10000PW[3,] <- scores32[3,]
    scores10000PW[4,] <- scores42[3,]
    scores10000PW[5,] <- scores43[3,]
    scores10000PW[6,] <- scores52[3,]

    # marginalize over number of observed and hidden variables
    score10000 <- array(0,dim=c(3,2))
    score10000[1,] <- colSums(scores10000DC)/n
    score10000[2,] <- colSums(scores10000AN)/n
    score10000[3,] <- colSums(scores10000PW)/n

    # --- write to file, before marginalization
    res <- list(s100DC=scores100DC,s1000DC=scores1000DC,s10000DC=scores10000DC,
                s100AN=scores100AN,s1000AN=scores1000AN,s10000AN=scores10000AN,
                s100PW=scores100PW,s1000PW=scores1000PW,s10000PW=scores10000PW)
    q_write <- q*100
    filename <- paste(path,"RESdisc",q_write,alg,nrep,"_",sig_write,".txt",
      sep="")
    writeListToFile(filename,res)

    return(list(s100=score100, s1000=score1000, s10000=score10000))

  } # end q=0.5

} # end function
