# This file contains functions to write data and knowledge to files, such that
# the TETRAD code can read it, as well as a function to read in the TETRAD-
# output
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



writeGraphToFile <- function(B, nobs, nhid, nrep, closure=NULL, graphFileName=
  "graph.txt", pathName="./") {

  # B: connection matrix
  # Example:
  # B = 1 0 0    x1(t+1) = x1(t)
  #     0 0 1    x2(t+1) = x3(t)
  #     1 1 0    x3(t+1) = x1(t) + x2(t)
  # nobs: number of observed variables (in 1 time slice)
  # nhid: number of hidden variables (in 1 time slice)
  # nrep: how often nodes will be repeated (= tau + 1 in paper)

  nobs_all <- nobs*nrep
  if (is.null(closure)) nhid_all <- nhid*nrep
  else nhid_all <- nhid*nrep + ncol(closure)

  ## get nodes string

  var_names <- NULL # for writing the edges

  str_var <- paste("X0",1,sep="")
  var_names <- paste("X0",1,sep="")

  for (i in 2:nobs_all) {
    if (i < 10) {
      str <- paste("X0",i,sep="")
    }
    else {
      str <- paste("X",i,sep="")
    }
    str_var <- paste(str_var,",",str,sep="")
    var_names <- c(var_names,str)
  }

  for (i in 1:nhid_all) {
    str <- paste("Latent(L",i,")",sep="")
    str_var <- paste(str_var,",",str,sep="")
    var_names <- c(var_names,paste("L",i,sep=""))
  }


  ## get edges string
  # small problem, have variables   4  1   but want to have  1  4
  #                                 5  2                     2  5
  #                                 6  3                     3  6
  # get it right with small trick

  B_all <- rep_B_sort(t(B),nrep-1,nobs=nobs)
  B_all[B_all==1] <- 2

  if (!is.null(closure)) {
    if (!ncol(closure)==0) {
      B_all <- add_closure_new(B_all,closure,nobs+nhid,nobs,nrep)
    }
  }

  for (i in 1:nrow(B_all)) {
    for (j in 1:ncol(B_all)) {
      if (B_all[i,j] == 1) B_all[j,i] <- 2
      else if (B_all[i,j] == 2) B_all[j,i] <- 1
    }
  }

  B_all[B_all==2] <- 0

  str_edge <- NULL

  # edges j -> i
  for (i in 1:(nobs_all+nhid_all)) {
    for (j in 1:(nobs_all+nhid_all)) {
      if (B_all[i,j]!=0) {
        str <- paste(var_names[j], "-->", var_names[i], sep="")
        if (is.null(str_edge)) str_edge <- str
        else str_edge <- paste(str_edge,",",str,sep="")
      }
    }
  }

  ## and finally write the whole string to the file
  str_final <- paste(str_var,",",str_edge,sep="")
  path <- paste(pathName,graphFileName,sep="")
  write(str_final,file=path,append=FALSE,sep="")

} # end function


writeDataToFile <- function(X, nrep, dataFileName="data.txt", pathName="./") {

  # X: data matrix with variables in columns, observations in rows
  # nrep: how often variables will be repeated (= tau + 1 in paper)

  nobs <- ncol(X)
  nsample <- nrow(X)
  nobs_all <- nobs*nrep

  neff <- nsample - nrep + 1
  cat("Number of samples:", neff, "\n")

  Y <- array(0,dim=c(neff,nobs_all))

  for (i in 1:neff) {
    temp <- X[i:(i+nrep-1),] # X1 earliest variable
    Y[i,] <- t(temp)
  }

  var_names <- array(0,dim=c(1,nobs_all))
  for (i in 1:nobs_all) {
    if (i < 10) {
      var_names[1,i] <- paste('X', 0, i, sep="")
    }
    else {
      var_names[1,i] <- paste('X', i, sep="")
    }
  }

  ## and finally write the whole string to the file
  path <- paste(pathName,dataFileName,sep="")
  write.table(var_names,file=path,append=FALSE,row.names=FALSE,col.names=FALSE)

  write.table(Y, file=path, append=TRUE, row.names=FALSE, col.names=FALSE)

} # end function


writeCovarianceToFile <- function(covMat, nsample, covFileName="cov.txt",
  pathName="./") {

  nobs_all <- dim(covMat)[1]

  path <- paste(pathName,covFileName,sep="")

  write("/covariance", file = path, append = FALSE)

  # sample size
  write(sprintf("%d",nsample), file = path, append = TRUE)

  # variable names
  var_names <- "X01"
  for (i in 2:nobs_all) {
    if (i < 10) {
      var_names <- paste(var_names, paste('X', 0, i, sep=""), sep=" ")
    }
    else {
      var_names <-paste(var_names,  paste('X', i, sep=""), sep=" ")
    }
  }

  write(var_names, file = path, append = TRUE)

  for (i in 1:nobs_all) {
    temp <- format(covMat[i,1],digits=20)
    for (j in index(2,i)) {
      temp <- paste(temp, format(covMat[i,j],digits=20), sep=" ")
    }
    write(temp, file = path, append = TRUE)
  }

}


writeKnowledgeToFile <- function(nobs,nrep,knowledgeFileName=
  "knowledge.txt",pathName="./") {

  # B: connection matrix
  # nobs: number of observed variables (in 1 time slice)
  # nrep: how often nodes will be repeated (= tau + 1 in paper)

  nobs_all <- nobs*nrep

  ## write tiers to file

  path <- paste(pathName,knowledgeFileName,sep="")
  if (nobs_all<10) {
    str <- paste(nrep-1," X", 0, nobs_all,sep="")
  }
  else {
    str <- paste(nrep-1," X",nobs_all,sep="")
  }
  write(str,file=path,append=FALSE,sep="")

  for (j in (nobs-1):1) {
    if (j+(nrep-1)*nobs < 10) {
      str <- paste(nrep-1," X",0,j+(nrep-1)*nobs,sep="")
    }
    else {
      str <- paste(nrep-1," X",j+(nrep-1)*nobs,sep="")
    }
    write(str,file=path,append=TRUE,sep="")
  }
  for (i in nrep:2) {
    for (j in nobs:1) {
      if (j+(i-2)*nobs < 10) {
        str <- paste(i-2," X",0,j+(i-2)*nobs,sep="")
      }
      else {
        str <- paste(i-2," X",j+(i-2)*nobs,sep="")
      }
      write(str,file=path,append=TRUE,sep="")
    }
  }
  
} # end function


callTetradUniversal <- function(type, pathName="./",
  tetradJarName="tetradcmd-4.3.9-18.jar", algName="fci",
  graphFileName="graph.txt", covFileName="cov.txt",
  dataFileName="data.txt", dataType="continuous", sig=0.01,
  knowledgeFileName="knowledge.txt", outFileName="resultFCI.txt",
  inclInstEffect=FALSE) {

  # type = "graph" (normal/standard DAG with hidden variables - or a tsgraph,
  #                but treated as normal graph, ie. no background knowledge)
  #        "data" (data generated from a normal/standard DAG with hidden vars -
  #                or from a tsgraph but treated as normal data)
  #        "tsgraph" (DAG model with HVs and repeating structure (over time))
  #        "tsdata" (data generated from a tsgraph)
  #        "graphbg" (DAG model with tiers as background knowledge)
  #        "databg" (date from DAG model with tiers as background knowledge)
  #        "cov"  (if continuous data are given in covariance form)
  # alg = "fci" (calls standard FCI)
  #       "cfci" (calls standard CFCI)
  #       "tsfci" (calls tsFCI)
  #       "tscfci" (calls tsCFCI)
  # datatype = "continuous"
  #            "discrete"
  # inclIE ... only needed for tsdata, otherwise ignored
  #     if FALSE then no instantaneous effects assumed, so inst. edges oriented
  #     as double headed arrows

  if (type == "graph") {
    # calls standard FCI/CFCI for graphs
    progcall <- paste("java -jar ", paste(pathName,tetradJarName,sep=""),
      " -graph ", paste(pathName,graphFileName,sep=""), " -algorithm ",
      algName, " -outfile ", paste(pathName,outFileName,sep=""), sep="")
  }

  if (type == "tsgraph") {
    # calls tsFCI/tsCFCI for graphs
    progcall <- paste("java -jar ", paste(pathName,tetradJarName,sep=""),
      " -graph ", paste(pathName,graphFileName,sep=""), " -algorithm ",
      algName, " -knowledge ", paste(pathName,knowledgeFileName,sep=""),
      " -inclInstEffect ", inclInstEffect, " -outfile ",
      paste(pathName,outFileName,sep=""), sep="")
  }

  if (type == "graphbg") {
    # calls standard FCI/CFCI for graphs including tiers in backgroundknowledge
    progcall <- paste("java -jar ", paste(pathName,tetradJarName,sep=""),
      " -graph ", paste(pathName,graphFileName,sep=""), " -algorithm ",
      algName, " -knowledge ", paste(pathName,knowledgeFileName,sep=""),
      " -outfile ", paste(pathName,outFileName,sep=""), sep="")
  }

  if (type == "data") {
    # calls standard FCI/CFCI for data
    progcall <- paste("java -jar ", paste(pathName,tetradJarName,sep=""),
      " -data ", paste(pathName,dataFileName,sep=""), " -datatype ", dataType,
      " -significance ", sig, " -algorithm ", algName, " -outfile ",
      paste(pathName,outFileName,sep=""), sep="")
  }

  if (type == "tsdata") {
    # calls tsFCI/tsCFCI for data
    progcall <- paste("java -jar ", paste(pathName,tetradJarName,sep=""),
      " -data ", paste(pathName,dataFileName,sep=""), " -datatype ", dataType,
      " -significance ", sig, " -algorithm ", algName, " -knowledge ",
      paste(pathName,knowledgeFileName,sep=""), " -inclInstEffect ",
      inclInstEffect, " -outfile ", paste(pathName,outFileName,sep=""),
      sep="")
  }

  if (type == "databg") {
    # calls standard FCI/CFCI for data including tiers in backgroundknowledge
    progcall <- paste("java -jar ", paste(pathName,tetradJarName,sep=""),
      " -data ", paste(pathName,dataFileName,sep=""), " -datatype ", dataType,
      " -significance ", sig, " -algorithm ", algName, " -knowledge ",
      paste(pathName,knowledgeFileName,sep=""), " -outfile ", 
      paste(pathName,outFileName,sep=""),  sep="")
  }

  if (type == "cov") {
    # calls (C)FCI or ts(C)FCI with covariance matrix as input
    if (algName=="fci" | algName=="cfci") {
      progcall <- paste("java -jar ", paste(pathName,tetradJarName,sep=""),
        " -covariance ", paste(pathName,covFileName,sep=""),
        " -significance ", sig, " -algorithm ", algName, " -outfile ",
        paste(pathName,outFileName,sep=""), sep="")
    }

    if (algName=="tsfci" | algName=="tscfci") {
      progcall <- paste("java -jar ", paste(pathName,tetradJarName,sep=""),
        " -covariance ", paste(pathName,covFileName,sep=""),
        " -significance ", sig, " -algorithm ", algName, " -knowledge ",
        paste(pathName,knowledgeFileName,sep=""), " -inclInstEffect ",
        inclInstEffect, " -outfile ", paste(pathName,outFileName,sep=""),
        sep="")
    }
  }

  system(progcall)

}


readFCIoutput <- function(pathName="./", readFileName="resultFCI.txt") {

  con <- file(paste(pathName,readFileName,sep=""))
  open(con)
  numNodes <- as.integer(readLines(con, n=1))
  numEdges <- as.integer(readLines(con, n=1))
  Nodes <- readLines(con, n=1) # do i actually need this?
  Edges <- readLines(con, n=1)
  close(con)

  Edges <- substr(Edges,2,nchar(Edges)-1) # delete first and last symbol
  Elist <- strsplit(Edges,",") # list of all edges

  pag <- array(0,dim=c(numNodes,numNodes))

  for (i in index(1,numEdges)) {
    temp <- Elist[[1]][i]
    if (i!=1) temp <- substr(temp,2,nchar(temp))
    temp <- strsplit(temp," ")

    node1 <- as.integer(substr(temp[[1]][1],2,nchar(temp[[1]][1])))
    node2 <- as.integer(substr(temp[[1]][3],2,nchar(temp[[1]][3])))
    
    mark1 <- substr(temp[[1]][2],1,1)
    mark2 <- substr(temp[[1]][2],3,3)

    if (mark1 == "-") pag[node2,node1] <- 1
    else if (mark1 == "<") pag[node2,node1] <- 2
    else if (mark1 == "o") pag[node2,node1] <- 3

    if (mark2 == "-") pag[node1,node2] <- 1
    else if (mark2 == ">") pag[node1,node2] <- 2
    else if (mark2 == "o") pag[node1,node2] <- 3

  }

  pag
  
} # end function