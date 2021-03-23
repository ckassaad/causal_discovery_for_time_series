# This file contains functions to run the simulations in "On Causal discovery
# from Time-Series Data using tsFCI", which is called in Simulations_graph.R,
# Simulations_data_cont.R and Simulations_data_disc.R, as well as functions to
# calculate the scores of an output of these three functions.
# This function calls TETRAD.
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


doSimulations_universal <- function(nobs, nhid, edgeperc, alg, nrep=4,
      nsim=100, dataType=NULL, nsample=NULL, sig=NULL, makeplot=FALSE,
      theseed=-1,bg=FALSE) {

  # simulates (ts)(C)FCI on discrete/continuous data or on graphs

  # --- needed for both data and graph
  # nobs ... number of observed variables
  # nhid ... number of hidden variables
  # edgeperc ... percentage with which each edge will be included, in [0,1]
  # alg ... string: "fci", "cfci", "tsfci", "tscfci"
  # nrep ... how often each node is included/repeated in the model (= tau+1 in paper)
  # nsim ... for how many models with the above parameters is the algorithm 
  #          called

  # --- only needed for data
  # nsample ... sample size
  # datatype ... string: "continuous" or "discrete"
  # sig ... significance value for independence tests

  # --- technical parameters, for both data and graph
  # makeplot ... if generating model and resulting pag should be printed
  # theseed ... which seed is used
  # bg ... if TRUE, then for (C)FCI background knowledge is included
  #        for ts(C)FCI this is ignored, since background knowledge is
  #        anyway included!!

  if (theseed == -1) {
    theseed <- sample(10000,1)
    set.seed(theseed)
  }
  else {
    set.seed(theseed)
  }

  cat("Selected seed was", theseed, ".\n")
#  wait(2)

  # seeds for separate runs - handy when using different sample sizes to get
  # the same generating models.
  allSeeds <- sample(nsim*100,nsim)

  res <- vector(mode="list",length=nsim)
  nvar <- nobs+nhid

  if (!is.null(dataType)) { # run on data
  
    if (is.null(sig)) return("Significance level needed.\n")
    if (is.null(nsample)) return("Sample size needed.\n")

    if (alg == "fci" | alg == "cfci") {
      if(!bg) temp_type <- "data"
      else temp_type <- "databg"
    }
    if (alg == "tsfci" | alg == "tscfci") temp_type <- "tsdata"
    if (nsample > 10000 & dataType == "continuous") temp_type <- "cov"
    if (nsample > 10000 & dataType == "discrete") {
      return("This is not feasible!")
    }

    print(temp_type)
  
    # generate data and call algorithm nsim times

    nburnin <- 1000 # how many samples are cut of from the beginnig

    if (dataType == "continuous") {
      for (i in 1:nsim) {

        cat("----------- This is run", i, ".---------------\n")

        set.seed(allSeeds[i])

        # matrix to graph in Figure 1
        # call with
        # thisseed <- 1 # to get result in Figure 2(c) use any of those: 1, 23, 41, 48, 53, 75
        # doSimulations_universal(2,1,0.25,"cfci",3,1,"continuous",100,0.01,TRUE,thisseed)
#        M1 <- matrix(0,3,3)
#        M1[1,1] <- 0.75
#        M1[2,c(1,3)] <- c(0.8,0.7)
#        M1[3,3] <- 0.88
#        temp <- genData_continuous(nobs,nhid,nsample,nburnin,edgeperc,M1)

        temp <- genData_continuous(nobs,nhid,nsample,nburnin,edgeperc)
        Xobs <- temp$X
        M1 <- temp$M1

        # define the connection matrix, for plot
        B <- M1
        B[B!=0] <- 1
        
        if (makeplot) temp_plot_fct(nobs, nhid, nrep, B)

        # call algorithm
        pag <- main_Tetrad_data(temp_type, t(Xobs), alg, nobs, nrep, sig,
          dataType)

        if (makeplot) plot_ts_pag(pag,nobs)

        res[[i]] <- list(M1=M1, pag=pag)

      }
    }

    if (dataType == "discrete") {
      for (i in 1:nsim) {

        cat("----------- This is run", i, ".---------------\n")

        set.seed(allSeeds[i])

        temp <- genData_discrete(nobs,nhid,nsample,nburnin,edgeperc)
        Xobs <- temp$X
        B <- temp$B
        CPT <- temp$CPT

        if (makeplot) temp_plot_fct(nobs, nhid, nrep, B)

        pag <- main_Tetrad_data(temp_type, t(Xobs), alg, nobs, nrep, sig, 
          dataType)

        if (makeplot) plot_ts_pag(pag,nobs)

        res[[i]] <- list(B=B, CPT=CPT, pag=pag)
      }
    }

    # write result to file
    ep_write <- edgeperc*100
    sig_write <- sig*100
    if (dataType=="continuous") dt_write <- "cont"
    else dt_write <- "disc"

    path <- "./simulation_results/"
    filename <- paste("DATA", dt_write, "_", alg, "_nobs", nobs, "nhid", nhid, 
        "edgep", ep_write, "nrep", nrep, "nsamp", nsample, "sig", sig_write,
        ".txt", sep="")
    writeListToFile(paste(path,filename, sep=""),res)

    # write generating models in table to load from matlab
    # generating models will be the same when only varying the sample since,
    # only write once
    if (nsample==100 && alg=="tscfci" && nrep==4 && sig==0.01) {
     filename <- paste("DATA", dt_write, "_", alg, "_nobs", nobs, "nhid", nhid, 
        "edgep", ep_write, ".txt", sep="")
      makeTableFromList(res, nobs, nhid, filename, path=path)
    }

  } 
#
  else { # run on graph

    if (alg == "fci" | alg == "cfci") {
      if(!bg) temp_type <- "graph"
      else temp_type <- "graphbg"
    }
    if (alg == "tsfci" | alg == "tscfci") temp_type <- "tsgraph"

    for (i in 1:nsim) {

      cat("----------- This is run", i, ".---------------\n")
  
      # generate randomly a connection matrix B - select parents
      # each edge is there with probability edgeperc

      set.seed(allSeeds[i])

      whichedge <- which(runif(nvar^2)<edgeperc)
      B <- array(0, dim=c(nvar,nvar))
      if (length(whichedge)!=0) B[whichedge] <- 1
#      B[1,c(2,3)] <- 1
#      B[2,c(2,3)] <- 1
#      B[3,c(1,3,4)] <- 1
#      B[4,c(3,4)] <- 1

      if (makeplot) temp_plot_fct(nobs, nhid, nrep, B)

      pag <- main_Tetrad_graph(temp_type, B, alg, nobs, nhid, nrep)

      if (makeplot) plot_ts_pag(pag,nobs)

      res[[i]] <- list(B=B,pag=pag)

    }

    # write result to file
    ep_write <- edgeperc*100

    path <- "./simulation_results/"
    filename <- paste("GRAPH_", alg, "_nobs", nobs, "nhid", nhid, "edgep",
        ep_write, "nrep", nrep, ".Rdata", sep="")
    writeListToFile(paste(path,filename, sep=""),res)

  }

  res

}




main_Tetrad_data <- function(type, Xobs, alg, nobs, nrep, sig, dataType, 
  inclIE=FALSE) {

  # runs (ts)(C)FCI for a given discrete/continuous data set

  # type ... what kind of input (time series or not, data or graph - see in the
  #          file "Tetrad_R_interact")
  #          "data", "tsdata", "graph", "tsgraph", "cov"
  # Xobs ... data matrix with variables in columns and observations in rows
  # alg ... string: "fci", "cfci", "tsfci", "tscfci"
  # nobs ... number of observed variables (columns of Xobs)
  # nrep ... how often each node is included/repeated in the model
  # sig ... significance value for independence tests
  # datatype ... string: "continuous" or "discrete"
  # inclIE ... only needed for tsdata, otherwise ignored
  #     if FALSE then no instantaneous effects assumed, so inst. edges oriented
  #     as double headed arrows

  if (dim(Xobs)[1] <= 10000) {
    # write data to file
    print("Write data to file.")
    writeDataToFile(Xobs, nrep)
  }
  else { # more than 10.000 samples (at least up to 1.000.000 can create all data first)
    if (dataType=="continuous") {
      # calculate covariance matrix
      n <- dim(Xobs)[1]
      covMat <- matrix(0,nobs*nrep,nobs*nrep)
      for (i in 1:(nobs*nrep)) {
        offseti <- floor((i-1)/nobs)
        indi <- i%%nobs
        if (indi==0) indi <- nobs
        tempi <- Xobs[(offseti+1):(n-nrep+offseti+1),indi,drop=FALSE]
        tempi <- tempi - mean(tempi)
        for (j in 1:i) {
          offsetj <- floor((j-1)/nobs)
          indj <- j%%nobs
          if (indj==0) indj <- nobs
          tempj <- Xobs[(offsetj+1):(n-nrep+offsetj+1),indj,drop=FALSE]
          tempj <- tempj-mean(tempj)
#          covMat[i,j] <- 1/(length(tempi)-1) * (t(tempi) %*% tempj)
          covMat[i,j] <- cov(tempi,tempj)
        }
      }
      # write covariance matrix to file
      writeCovarianceToFile(covMat, n)   
    }
    else {
      return("Discrete data with more than 10.000 data not feasible.")
    }
  }

  # write background knowledge to file
  if (type == "tsdata" || type == "databg" ||
     (type == "cov" & (alg == "tsfci" | alg == "tscfci")) ) {
    print("Write knowledge to file.")
    writeKnowledgeToFile(nobs, nrep)
  }

  # run alg from Tetrad, save output PAG to "resultFCI.txt"
  print("Call Tetrad.")

  callTetradUniversal(type, algName=alg, dataFileName="data.txt", dataType=dataType, sig=sig,
    inclInstEffect=inclIE)

  # read in results from Fci
  print("Read in (C)FCI output.")
  pag <- readFCIoutput()

  pag

}


main_Tetrad_graph <- function(type, B, alg, nobs, nhid, nrep,inclIE=FALSE) {

  # runs (ts)(C)FCI for a given graph (generating model with perfect oracle)

  # type ... what kind of input (time series or not, data or graph)
  #          "data", "tsdata", "graph", "tsgraph"
  # B ... connection matrix of graph
  # alg ... string: "fci", "cfci", "tsfci", "tscfci"
  # nobs ... number of observed variables
  # nhid ... number of hidden variables
  # nrep ... how often each node is included/repeated in the model
  # inclIE ... only needed for tsdata, otherwise ignored
  #     if FALSE then no instantaneous effects assumed, so inst. edges oriented
  #     as double headed arrows

  closure <- get_closure(B) # common ancestors of earliest variables

  # write graph to file
  print("write graph to file")
  writeGraphToFile(B,nobs,nhid,nrep,closure)

  # write background knowledge to file
  if (type == "tsgraph" || type == "databg") {
    print("Write knowledge to file.")
    writeKnowledgeToFile(nobs, nrep)
  }

  # run alg from Tetrad, save output PAG to "resultFCI.txt" in folder rfci
  print("Call Tetrad.")

  callTetradUniversal(type, algName=alg, inclInstEffect=inclIE)

  # read in results from Fci
  print("Read in (C)FCI output.")
  pag <- readFCIoutput()

  pag

}


temp_plot_fct <- function(nobs, nhid, nrep, B) {

    # plots the time series of the generating model

    nvar <- nobs+nhid

    B_inst <- array(0,dim=c(nvar,nvar))
    ngvec <- c(rep(FALSE,nobs),rep(TRUE,nhid))
    nodelabels <- NULL
    for (i in 1:nrep) {
      temp <- c(1:nobs+(i-1)*nobs,nobs*nrep+index(1,nhid)+(i-1)*nhid)
      nodelabels <- c(nodelabels,temp)
    }

    closure <- get_closure(B) # common ancestors of earliest variables
    if (ncol(closure)==0){
      plot_timeseries(list(B_inst,B),ngvec=ngvec,nrepeat=nrep-2,nodelabels=
        nodelabels)
    }
    else {
      plot_timeseries_closure(list(B_inst,B),closure,ngvec=ngvec,nrepeat=
        nrep-2,hidden=TRUE,nobs=nobs)
    }

}


checkSimulationOutput_scores <- function(simOut,nobs,nrep=4) {

  # calculates all three scores.
  # simOut is (a list) of outputs of a simulation of (ts)(C)FCI

  nsim <- length(simOut)

  scores1 <- array(0,dim=c(nsim,2)) # direct cause score
  scores2 <- array(0,dim=c(nsim,2)) # ancestor score
  scores3 <- array(0,dim=c(nsim,2)) # pairwise score

  for (i in 1:nsim) {

    if (is.null(simOut[[i]]$B)) {
      M1 <- simOut[[i]]$M1
      B <- M1
      B[B!=0] <- 1
    }
    else B <- simOut[[i]]$B

    pag <- simOut[[i]]$pag

    postproc <- TRUE
    if (postproc) {
      # in (C)FCI orient all arrows forward in time as postprocessing
      # in ts(C)FCI this is anyway given in the output
      temp <- pag
      temp[lower.tri(temp,TRUE)] <- 0
      temp[temp==3] <- 2
      pag[upper.tri(pag)] <- temp[upper.tri(temp)]
    }

    makeplot <- FALSE # makeplot <- TRUE
    if (makeplot) temp_plot_fct(nobs, dim(B)[1]-nobs, nrep, B)
    if (makeplot) plot_ts_pag(pag,nobs)

    scores1[i,] <- compareGenModelToPAG_direct(B,pag,nrep,nobs)$scores
    scores2[i,] <- compareGenModelToPAG_ancestor(B,pag,nrep,nobs)$scores
    scores3[i,] <- compareGenModelToPag_pairwise(B,pag,nrep,nobs)$scores

  }

#  print(scores1)
#  print(scores2)  
#  print(scores3)

  score_sum <- array(0,dim=c(3,2))
  score_sum[1,1] <- sum(scores1[,1])/(nsim-sum(is.na(scores1[,2])))
  score_sum[1,2] <- sum(scores1[,2],na.rm=TRUE)/(nsim-sum(is.na(scores1[,2])))
  score_sum[2,1] <- sum(scores2[,1])/(nsim-sum(is.na(scores2[,2])))
  score_sum[2,2] <- sum(scores2[,2],na.rm=TRUE)/(nsim-sum(is.na(scores2[,2])))
  score_sum[3,1] <- sum(scores3[,1])/(nsim-sum(is.na(scores3[,2])))
  score_sum[3,2] <- sum(scores3[,2],na.rm=TRUE)/(nsim-sum(is.na(scores3[,2])))

  list(scores=score_sum)

}



GrangerSim <- function(simOutGraph,nobs,nrep=4) {

  # SimOut ... Output of a Simulation of tsFCI using graphs (see above)
  # looks for each graph in SimOut which variables Granger-cause each other in
  # infinite sample size

  nsim <- length(simOutGraph) # number of different graphs

  scores1 <- array(0,dim=c(nsim,2)) # direct
  scores2 <- array(0,dim=c(nsim,2)) # ancestor
  scores3 <- array(0,dim=c(nsim,2)) # pairwise

  for (i in 1:nsim) {

    temp <- Granger_inf(simOutGraph[[i]]$B,nrep,nobs)$scores
    scores1[i,] <- temp[1,]
    scores2[i,] <- temp[2,]
    scores3[i,] <- temp[3,]

  }

  score_sum <- array(0,dim=c(3,2))
  score_sum[1,1] <- sum(scores1[,1])/(nsim-sum(is.na(scores1[,2])))
  score_sum[1,2] <- sum(scores1[,2],na.rm=TRUE)/(nsim-sum(is.na(scores1[,2])))
  score_sum[2,1] <- sum(scores2[,1])/(nsim-sum(is.na(scores2[,2])))
  score_sum[2,2] <- sum(scores2[,2],na.rm=TRUE)/(nsim-sum(is.na(scores2[,2])))
  score_sum[3,1] <- sum(scores3[,1])/(nsim-sum(is.na(scores3[,2])))
  score_sum[3,2] <- sum(scores3[,2],na.rm=TRUE)/(nsim-sum(is.na(scores3[,2])))

  list(scores=score_sum)

}


Granger_inf <- function(B,nrep,nobs) {

  nvar <- dim(B)[1]

  B1 <- rep_B_sort(t(B),nrep-1,nobs=nobs)
  B1[B1==1] <- 2

  closure <- get_closure(B)
  B1 <- add_closure_new(B1,closure,nvar,nobs,nrep)

  for (i in 1:nrow(B1)) {
    for (j in 1:ncol(B1)) {
      if (B1[i,j] == 1) B1[j,i] <- 2
      else if (B1[i,j] == 2) B1[j,i] <- 1
    }
  }

  grangerGraph <- array(0,dim=c(nobs*nrep,nobs*nrep))

  # go through all the pairs with one variable in latest time layer
  for (x in (nobs*(nrep-1)+1):(nobs*nrep)) { # later in time 

    MB <- findMarkovBlanket(x,B1,nobs,nvar,nrep)
   
    for (y in MB) { # earlier in time
      if ( y <= (nobs*nrep)) { # y observed
        grangerGraph[x,y] <- 1
        grangerGraph[y,x] <- 2
        
        # fill in for all similar pairs
        for (i in 1:(nrep-2)) {
          aa <- i*nobs
          if (y-aa <= 0) break
          grangerGraph[x-aa,y-aa] <- 1
          grangerGraph[y-aa,x-aa] <- 2
        }
      }
    }

  }

  makeplot <- FALSE
  if (makeplot) {
    print("Plot Granger Graph.")
    plot_ts_pag(grangerGraph,nobs)
  }

  score <- array(0,dim=c(3,2))

  score[1,] <- compareGenModelToGranger_direct(B,grangerGraph,nrep,nobs)$score
  score[2,] <- compareGenModelToGranger_ancestor(B,grangerGraph,nrep,nobs)$score
  score[3,] <- compareGenModelToGranger_pairwise(B,grangerGraph,nrep,nobs)$score

  return(list(GG=grangerGraph,scores=score))

}



writeListToFile <- function(filename,temp) {
  con <- file(filename)
  open(con,"w")
  serialize(temp,con)
  close(con)
}

readListFromFile <- function(filename){
  con <- file(filename,"r")
  simOut <- unserialize(con)
  close(con)
  simOut
}


makeTableFromList <- function(simout, nobs, nhid, filename, path=
  "./simulation_results/") {

  # read out file form simulations and save matrices in a table -> to read in
  # from matlab
  # path = "./simulation_results/"
  # filename = "DATAcont_tscfci_nobs2nhid1edgep25nrep4nsamp100sig1.txt"
  # first row of table contains number of observed and hidden variables

  dims <- dim(simout[[1]][[1]])[1]
  
  conmat <- array(0,dim=c(length(simout)*dims+1,dims))
  conmat[1,1:2] <- c(nobs,nhid)
  k <- 2
  
  for (i in 1:length(simout)) {
    conmat[k:(k+dims-1),] <- simout[[i]][[1]]
    k <- k+dims
  }

  writeTo <- paste(path, "Table_", filename, sep="")

  write.table(conmat, file=writeTo, row.names=FALSE, col.names=FALSE)

}


