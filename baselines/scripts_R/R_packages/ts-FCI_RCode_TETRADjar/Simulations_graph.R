# This file contains the commands to run the simulations in the infinite sample
# size limit in "On Causal discovery from Time-Series Data using tsFCI"
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
# -------------------- simulations in infinite sample size ------------------- #
# --------------------------------------------------------------------------- #

simInfSamp <- function(q=0.25, nrep=4, nsim=100) {

  # nrep = tau+1 in paper

  path <- "./simulation_results/" # where results will be saved

  if (q!=0.25 & q!=0.5) {
    return("Input parameter q must either be 0.25 or 0.5.")
  }

# ------------------------- edgepercentage q = 0.25 ------------------------- #
  if (q==0.25) {

    # this does "nsim" simulations for each setting of parameters for tsfci,
    # granger and fci in infinite sample size limit
    # saves generating models and resulting pags for each setting in a file in 
    #   the folder simulation_results (for tsfci and fci only)
    # also saves a summary of the socres in the folder simulation_results
    # can load files using readListFromFile(filenname) - see main_tetrad_fci.R

    seeds <- c(3636,8309,9035,1502,2029,908)
    n <- length(seeds) # number of different seetings for nobs and nhid

    # --- tsfci

    simOut21 <- doSimulations_universal(2,1,0.25,"tsfci",nrep,nsim,theseed=seeds[1])
    simOut22 <- doSimulations_universal(2,2,0.25,"tsfci",nrep,nsim,theseed=seeds[2])
    simOut31 <- doSimulations_universal(3,1,0.25,"tsfci",nrep,nsim,theseed=seeds[3])
    simOut32 <- doSimulations_universal(3,2,0.25,"tsfci",nrep,nsim,theseed=seeds[4])
    simOut42 <- doSimulations_universal(4,2,0.25,"tsfci",nrep,nsim,theseed=seeds[5])
    simOut33 <- doSimulations_universal(3,3,0.25,"tsfci",nrep,nsim,theseed=seeds[6])

    scores21 <- checkSimulationOutput_scores(simOut21,2,nrep)$scores
    scores22 <- checkSimulationOutput_scores(simOut22,2,nrep)$scores
    scores31 <- checkSimulationOutput_scores(simOut31,3,nrep)$scores
    scores32 <- checkSimulationOutput_scores(simOut32,3,nrep)$scores
    scores42 <- checkSimulationOutput_scores(simOut42,4,nrep)$scores
    scores33 <- checkSimulationOutput_scores(simOut33,3,nrep)$scores

    # direct-cause scores
    scoresTsfciDC <- array(0,dim=c(n,2))
    scoresTsfciDC[1,] <- scores21[1,]
    scoresTsfciDC[2,] <- scores22[1,]
    scoresTsfciDC[3,] <- scores31[1,]
    scoresTsfciDC[4,] <- scores32[1,]
    scoresTsfciDC[5,] <- scores42[1,]
    scoresTsfciDC[6,] <- scores33[1,]

    # ancestor scores
    scoresTsfciAN <- array(0,dim=c(n,2))
    scoresTsfciAN[1,] <- scores21[2,]
    scoresTsfciAN[2,] <- scores22[2,]
    scoresTsfciAN[3,] <- scores31[2,]
    scoresTsfciAN[4,] <- scores32[2,]
    scoresTsfciAN[5,] <- scores42[2,]
    scoresTsfciAN[6,] <- scores33[2,]

    # pairwise scores
    scoresTsfciPW <- array(0,dim=c(n,2))
    scoresTsfciPW[1,] <- scores21[3,]
    scoresTsfciPW[2,] <- scores22[3,]
    scoresTsfciPW[3,] <- scores31[3,]
    scoresTsfciPW[4,] <- scores32[3,]
    scoresTsfciPW[5,] <- scores42[3,]
    scoresTsfciPW[6,] <- scores33[3,]

    # marginalize over number of observed and hidden variables
    scoreTsfci <- array(0,dim=c(3,2))
    scoreTsfci[1,] <- colSums(scoresTsfciDC)/n
    scoreTsfci[2,] <- colSums(scoresTsfciAN)/n
    scoreTsfci[3,] <- colSums(scoresTsfciPW)/n

    # --- Granger

    # check Granger causal relationships for same graphs as used in the simulations
    # of tsFCI

    scores21 <- GrangerSim(simOut21,2,nrep)$scores
    scores22 <- GrangerSim(simOut22,2,nrep)$scores
    scores31 <- GrangerSim(simOut31,3,nrep)$scores
    scores32 <- GrangerSim(simOut32,3,nrep)$scores
    scores42 <- GrangerSim(simOut42,4,nrep)$scores
    scores33 <- GrangerSim(simOut33,3,nrep)$scores

    # direct cause score
    scoresGrDC <- array(0,dim=c(n,2))
    scoresGrDC[1,] <- scores21[1,]
    scoresGrDC[2,] <- scores22[1,]
    scoresGrDC[3,] <- scores31[1,]
    scoresGrDC[4,] <- scores32[1,]
    scoresGrDC[5,] <- scores42[1,]
    scoresGrDC[6,] <- scores33[1,]

    # ancestor score
    scoresGrAN <- array(0,dim=c(n,2))
    scoresGrAN[1,] <- scores21[2,]
    scoresGrAN[2,] <- scores22[2,]
    scoresGrAN[3,] <- scores31[2,]
    scoresGrAN[4,] <- scores32[2,]
    scoresGrAN[5,] <- scores42[2,]
    scoresGrAN[6,] <- scores33[2,]

    # pairwise score
    scoresGrPW <- array(0,dim=c(n,2))
    scoresGrPW[1,] <- scores21[3,]
    scoresGrPW[2,] <- scores22[3,]
    scoresGrPW[3,] <- scores31[3,]
    scoresGrPW[4,] <- scores32[3,]
    scoresGrPW[5,] <- scores42[3,]
    scoresGrPW[6,] <- scores33[3,]

    # marginalize over number of observed and hidden variables
    scoreGr <- array(0,dim=c(3,2))
    scoreGr[1,] <- colSums(scoresGrDC)/n
    scoreGr[2,] <- colSums(scoresGrAN)/n
    scoreGr[3,] <- colSums(scoresGrPW)/n

    # --- original FCI

    # run original FCI on same graphs as used in the simulations of tsFCI by 
    # setting the same seed (last number in input)
    simOut21 <- doSimulations_universal(2,1,0.25,"fci",nrep,nsim,theseed=seeds[1])
    simOut22 <- doSimulations_universal(2,2,0.25,"fci",nrep,nsim,theseed=seeds[2])
    simOut31 <- doSimulations_universal(3,1,0.25,"fci",nrep,nsim,theseed=seeds[3])
    simOut32 <- doSimulations_universal(3,2,0.25,"fci",nrep,nsim,theseed=seeds[4])
    simOut42 <- doSimulations_universal(4,2,0.25,"fci",nrep,nsim,theseed=seeds[5])
    simOut33 <- doSimulations_universal(3,3,0.25,"fci",nrep,nsim,theseed=seeds[6])

    scores21 <- checkSimulationOutput_scores(simOut21,2,nrep)$scores
    scores22 <- checkSimulationOutput_scores(simOut22,2,nrep)$scores
    scores31 <- checkSimulationOutput_scores(simOut31,3,nrep)$scores
    scores32 <- checkSimulationOutput_scores(simOut32,3,nrep)$scores
    scores42 <- checkSimulationOutput_scores(simOut42,4,nrep)$scores
    scores33 <- checkSimulationOutput_scores(simOut33,3,nrep)$scores

    # direct-cause scores
    scoresFciDC <- array(0,dim=c(n,2))
    scoresFciDC[1,] <- scores21[1,]
    scoresFciDC[2,] <- scores22[1,]
    scoresFciDC[3,] <- scores31[1,]
    scoresFciDC[4,] <- scores32[1,]
    scoresFciDC[5,] <- scores42[1,]
    scoresFciDC[6,] <- scores33[1,]

    # ancestor scores
    scoresFciAN <- array(0,dim=c(n,2))
    scoresFciAN[1,] <- scores21[2,]
    scoresFciAN[2,] <- scores22[2,]
    scoresFciAN[3,] <- scores31[2,]
    scoresFciAN[4,] <- scores32[2,]
    scoresFciAN[5,] <- scores42[2,]
    scoresFciAN[6,] <- scores33[2,]

    # pairwise scores
    scoresFciPW <- array(0,dim=c(n,2))
    scoresFciPW[1,] <- scores21[3,]
    scoresFciPW[2,] <- scores22[3,]
    scoresFciPW[3,] <- scores31[3,]
    scoresFciPW[4,] <- scores32[3,]
    scoresFciPW[5,] <- scores42[3,]
    scoresFciPW[6,] <- scores33[3,]

    # marginalize over number of observed and hidden variables
    scoreFci <- array(0,dim=c(3,2))
    scoreFci[1,] <- colSums(scoresFciDC)/n
    scoreFci[2,] <- colSums(scoresFciAN)/n
    scoreFci[3,] <- colSums(scoresFciPW)/n

    # --- write to file, before marginalization
    res <- list(sTsfciDC=scoresTsfciDC, sGrDC=scoresGrDC, sFciDC=scoresFciDC,
                sTsfciAN=scoresTsfciAN, sGrAN=scoresGrAN, sFciAN=scoresFciAN,
                sTsfciPW=scoresTsfciPW, sGrPW=scoresGrPW, sFciPW=scoresFciPW)
    q_write <- q*100
    filename <- paste(path,"RESinf",q_write,"_",nrep,".txt",sep="")
    writeListToFile(filename,res)

    return(list(sTsfci=scoreTsfci, sGranger=scoreGr, sFci=scoreFci))

  } # end q=0.25

# ------------------------- edgepercentage q = 0.50 ------------------------- #
  if (q==0.5) {

    # this does "nsim" simulations for each setting of parameters for tsfci,
    # granger and fci in infinite sample size limit
    # saves generating models and resulting pags for each setting in a file in 
    #   the folder simulation_results (for tsfci and fci only)
    # also saves a summary of the socres in the folder simulation_results
    # can load files using readListFromFile(filenname) - see main_tetrad_fci.R

    seeds <- c(9879,6789,7562,5833)
    n <- length(seeds) # number of different seetings for nobs and nhid

    # --- tsFCI

    simOut21 <- doSimulations_universal(2,1,0.5,"tsfci",nrep,nsim,theseed=seeds[1])
    simOut22 <- doSimulations_universal(2,2,0.5,"tsfci",nrep,nsim,theseed=seeds[2])
    simOut31 <- doSimulations_universal(3,1,0.5,"tsfci",nrep,nsim,theseed=seeds[3])
    simOut32 <- doSimulations_universal(3,2,0.5,"tsfci",nrep,nsim,theseed=seeds[4])

    scores21 <- checkSimulationOutput_scores(simOut21,2,nrep)$scores
    scores22 <- checkSimulationOutput_scores(simOut22,2,nrep)$scores
    scores31 <- checkSimulationOutput_scores(simOut31,3,nrep)$scores
    scores32 <- checkSimulationOutput_scores(simOut32,3,nrep)$scores

    # direct-cause scores
    scoresTsfciDC <- array(0,dim=c(n,2))
    scoresTsfciDC[1,] <- scores21[1,]
    scoresTsfciDC[2,] <- scores22[1,]
    scoresTsfciDC[3,] <- scores31[1,]
    scoresTsfciDC[4,] <- scores32[1,]

    # ancestor scores
    scoresTsfciAN <- array(0,dim=c(n,2))
    scoresTsfciAN[1,] <- scores21[2,]
    scoresTsfciAN[2,] <- scores22[2,]
    scoresTsfciAN[3,] <- scores31[2,]
    scoresTsfciAN[4,] <- scores32[2,]

    # pairwise scores
    scoresTsfciPW <- array(0,dim=c(n,2))
    scoresTsfciPW[1,] <- scores21[3,]
    scoresTsfciPW[2,] <- scores22[3,]
    scoresTsfciPW[3,] <- scores31[3,]
    scoresTsfciPW[4,] <- scores32[3,]

    # marginalize over number of observed and hidden variables
    scoreTsfci <- array(0,dim=c(3,2))
    scoreTsfci[1,] <- colSums(scoresTsfciDC)/n
    scoreTsfci[2,] <- colSums(scoresTsfciAN)/n
    scoreTsfci[3,] <- colSums(scoresTsfciPW)/n

    # --- Granger

    # check Granger causal relationships for same graphs as used in the simulations
    # of tsFCI

    scores21 <- GrangerSim(simOut21,2,nrep)$scores
    scores22 <- GrangerSim(simOut22,2,nrep)$scores
    scores31 <- GrangerSim(simOut31,3,nrep)$scores
    scores32 <- GrangerSim(simOut32,3,nrep)$scores

    # direct cause score
    scoresGrDC <- array(0,dim=c(n,2))
    scoresGrDC[1,] <- scores21[1,]
    scoresGrDC[2,] <- scores22[1,]
    scoresGrDC[3,] <- scores31[1,]
    scoresGrDC[4,] <- scores32[1,]

    # ancestor score
    scoresGrAN <- array(0,dim=c(n,2))
    scoresGrAN[1,] <- scores21[2,]
    scoresGrAN[2,] <- scores22[2,]
    scoresGrAN[3,] <- scores31[2,]
    scoresGrAN[4,] <- scores32[2,]

    # pairwise score
    scoresGrPW <- array(0,dim=c(n,2))
    scoresGrPW[1,] <- scores21[3,]
    scoresGrPW[2,] <- scores22[3,]
    scoresGrPW[3,] <- scores31[3,]
    scoresGrPW[4,] <- scores32[3,]

    # marginalize over number of observed and hidden variables
    scoreGr <- array(0,dim=c(3,2))
    scoreGr[1,] <- colSums(scoresGrDC)/n
    scoreGr[2,] <- colSums(scoresGrAN)/n
    scoreGr[3,] <- colSums(scoresGrPW)/n

    # --- original FCI

    # run original FCI on same graphs as used in the simulations of tsFCI by 
    # setting the same seed (last number in input)
    simOut21 <-doSimulations_universal(2,1,0.5,"fci",nrep,nsim,theseed=seeds[1])
    simOut22 <-doSimulations_universal(2,2,0.5,"fci",nrep,nsim,theseed=seeds[2])
    simOut31 <-doSimulations_universal(3,1,0.5,"fci",nrep,nsim,theseed=seeds[3])
    simOut32 <-doSimulations_universal(3,2,0.5,"fci",nrep,nsim,theseed=seeds[4])

    scores21 <- checkSimulationOutput_scores(simOut21,2,nrep)$scores
    scores22 <- checkSimulationOutput_scores(simOut22,2,nrep)$scores
    scores31 <- checkSimulationOutput_scores(simOut31,3,nrep)$scores
    scores32 <- checkSimulationOutput_scores(simOut32,3,nrep)$scores

    # direct-cause scores
    scoresFciDC <- array(0,dim=c(n,2))
    scoresFciDC[1,] <- scores21[1,]
    scoresFciDC[2,] <- scores22[1,]
    scoresFciDC[3,] <- scores31[1,]
    scoresFciDC[4,] <- scores32[1,]

    # ancestor scores
    scoresFciAN <- array(0,dim=c(n,2))
    scoresFciAN[1,] <- scores21[2,]
    scoresFciAN[2,] <- scores22[2,]
    scoresFciAN[3,] <- scores31[2,]
    scoresFciAN[4,] <- scores32[2,]

    # pairwise scores
    scoresFciPW <- array(0,dim=c(n,2))
    scoresFciPW[1,] <- scores21[3,]
    scoresFciPW[2,] <- scores22[3,]
    scoresFciPW[3,] <- scores31[3,]
    scoresFciPW[4,] <- scores32[3,]

    # marginalize over number of observed and hidden variables
    scoreFci <- array(0,dim=c(3,2))
    scoreFci[1,] <- colSums(scoresFciDC)/n
    scoreFci[2,] <- colSums(scoresFciAN)/n
    scoreFci[3,] <- colSums(scoresFciPW)/n

    # --- write to file, before marginalization
    res <- list(sTsfciDC=scoresTsfciDC, sGrDC=scoresGrDC, sFciDC=scoresFciDC,
                sTsfciAN=scoresTsfciAN, sGrAN=scoresGrAN, sFciAN=scoresFciAN,
                sTsfciPW=scoresTsfciPW, sGrPW=scoresGrPW, sFciPW=scoresFciPW)
    q_write <- q*100
    filename <- paste(path,"RESinf",q_write,"_",nrep,".txt",sep="")
    writeListToFile(filename,res)

    return(list(sTsfci=scoreTsfci, sGranger=scoreGr, sFci=scoreFci))

  } # end q=0.5

} # end function
