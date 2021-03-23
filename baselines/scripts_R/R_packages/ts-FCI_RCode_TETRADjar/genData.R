# This file contains functions to generate discrete and continuous data from
# a 1-lag model
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



genData_continuous <- function(nobs,nhid,nsample,nburnin,edgeperc,M1=NULL) {

  nvar <- nobs + nhid
  ntotal <- nsample + nburnin

  if (is.null(M1)) {
    # get stationary VAR-model; linear, Gaussian errors
    whichedge <- which(runif(nvar^2)<edgeperc) # decide which edges are present
    M1 <- genM1_random(nvar, nedge=length(whichedge), whichedge=whichedge)
    while (!all(abs(eigen(M1)$value)<1)) {
      M1 <- genM1_random(nvar, nedge=length(whichedge), whichedge=whichedge)
    }
  }

#  M1 <- array(0,dim=c(3,3))
#  M1[1,c(1,2,3)] <- c(0.72,0.97, -0.76)
#  M1[2,c(1,2,3)] <- c(-0.56,-0.58,1.22)
#  M1[3,c(1,3)] <- c(-0.89,0.97)
#  M1[4,c(1,3)] <- c(-0.62, 0.80)

#  M1[1,c(1)] <- c(0.72)
#  M1[2,c(3)] <- c(1.22)
#  M1[3,c(1)] <- c(-0.89)

#  M1 <- array(0,dim=c(3,3))
#  M1[1,1] <- 0.9
#  M1[2,3] <- 0.85
#  M1[3,1] <- 0.96

#  print(M1)
  cat("abs. value of EV's of M1", abs(eigen(M1)$value), "\n")

  c <- rep(0,nvar)

  # generate Data
  ntotal <- nsample+nburnin
  X <- array(0,dim=c(nvar,ntotal))
  err <- array(0,dim=c(nvar,ntotal))

  err[,1] <- rnorm(nvar,sd=1)
  X[,1] <- c + err[,1]
  for (i in 2:(ntotal)) {
    err[,i] <- rnorm(nvar,sd=1)
    X[,i] <- c + M1%*%X[,i-1] + err[,i]
  }

  # cut out beginning of the time series
  X <- X[,(nburnin+1):ntotal]
  # err <- err[,(nburnin+1):ntotal]

  X <- X[1:nobs,] # observed time series data

  list(X=X, M1=M1)

}


genData_discrete <- function(nobs,nhid,nsample,nburnin,edgeperc, B=NULL, 
  CPT=NULL) {

  nvar <- nobs + nhid

  if (is.null(B) & is.null(CPT)) {

    # generate randomly a connection matrix B - select parents
    whichedge <- which(runif(nvar^2)<edgeperc)

    B <- array(0, dim=c(nvar,nvar))
    if (length(whichedge)!=0) B[whichedge] <- 1
#    print(B)

    CPT <- vector(mode="list", length=nvar)

    # parents, defined by B
    pa <- vector(mode="list", length=nvar)
    for (x in 1:nvar) pa[[x]] <- which(B[x,]!=0)

    # choose number of discrete states per variable
    nstates <- vector(mode="integer",length=nvar)
    for (i in 1:nvar) {
      nstates[i] <- 2 # sample(2:4,1) #2 # use binary only
    }

#    print(nstates)

    # construct the CPTs
    for (x in 1:nvar) {

      dims <- nstates[x]
      temp <- 1
      for (y in pa[[x]]) temp <- temp*nstates[y]
      dims <- c(dims,temp)

      CPT[[x]] <- array(0,dim=dims)

      # this only works for binary variables
      if (length(pa[[x]])==0) {
#        temp <- runif(nstates[x],min=0.1,max=0.9)
#        temp <- temp/sum(temp) # distribution of x given 1 parents constell.
        temp1 <- runif(1)
        temp <- array(0,dim=c(1,2))
        if (temp1 < 0.5) {
          temp[1] <- runif(1,min=0.05, max=0.45)
          temp[2] <- 1-temp[1]
        }
        else {
          temp[1] <- runif(1,min=0.55, max=0.95)
          temp[2] <- 1-temp[1]
        }
        CPT[[x]] <- temp
        dim(CPT[[x]]) <- c(nstates[x],1)
      }
      else {
        all_dims <- dims[2]

        for (i in 1:all_dims) {
#          temp <- runif(nstates[x],min=0.1,max=0.9)
#          temp <- temp/sum(temp) # distribution of x given one parents constell.
          temp1 <- runif(1)
          temp <- array(0,dim=c(1,2))
          if (temp1 < 0.5) {
            temp[1] <- runif(1,min=0.05, max=0.45)
            temp[2] <- 1-temp[1]
          }
          else {
            temp[1] <- runif(1,min=0.55, max=0.95)
            temp[2] <- 1-temp[1]
          }
          CPT[[x]][,i] <- temp
        }
      }

    } # end construct CPTs

  }

  print(CPT)

  # generate data        
  ntotal <- nsample+nburnin
  X <- array(0,dim=c(nvar,ntotal))

  for (x in 1:nvar) X[x,1] <- sample(nstates[x],1)

  for (i in 2:(ntotal)) {
    for (x in 1:nvar) {
      pa_val <- X[pa[[x]],i-1]
      if (!(length(pa_val) == 0)) { # if node x has parents - get index of
                                    # of parents constellation
        pa_states <- nstates[pa[[x]]]
        # get the index of the parents constellation
        ind <- pa_val[1]
        if (length(pa_val)>1) {
          for (j in length(pa_val):2) {
            ind <- ind + (pa_val[j]-1)*prod(pa_states[(j-1):1])
          }
        }
      }
      else {
        ind <- 1
      }

      temp <- runif(1)
      x_val <- min(which(temp < cumsum(CPT[[x]][,ind]) ))
      X[x,i] <- x_val
    }
  }

  X <- X[,(nburnin+1):ntotal]
  X <- X[1:nobs,] # observed time series data

  list(X=X, B=B, CPT=CPT)

}



genM1_random <- function(dims, nedge=NULL, whichedge=NULL) {

  M1 <- array(0,dim=c(dims,dims))

  if (is.null(nedge)) {
    prob <- NULL
    # only relatively sparse graphs (total: dims^2)
    for (i in 1:(dims+3)) prob <- c(prob,choose(dims^2,i))
    cumprob <- cumsum(prob/sum(prob))
    nedge <- min(which(runif(1)<cumprob))
  }

  if (nedge==0) return(M1)

  if (is.null(whichedge)) {
    whichedge <- sample(dims^2,nedge)
  }

  for (i in 1:nedge) {
    r <- ceiling(runif(1)*2)
    M1[whichedge[i]] <- round(runif(1,min=0.2,max=1.5)*(-1)^r, digits=2)
  }

  M1

}
