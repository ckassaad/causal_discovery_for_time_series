# This file contains functions to calculate the direct-cause, ancestor and
# pairwise scores
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



compareGenModelToPAG_direct <- function(B,pag,nrep,nobs) {

  # "direct-cause" score
  # B ... the connection matrix of the generating model over all variables
  # pag ... the PAG returned by (ts)(C)FCI over ALL observed variales

  # for each order pair (y,x) with y earlier in time than x, and x in last time
  # slice, check what decision we got from the PAG and compare it to the
  # generating model
  # 1 ... y direct cause of x = causal effect
  # 2 ... y no direct cause of x = no causal effect
  # 3 ... don't know (circles)

  # for comparison to generating model: put connection matrix nrep times in
  # one big matrix, similar to pag, with 1 and 2 as endpoints
  B1 <- rep_B_sort(t(B),nrep-1,nobs=nobs)
  B1[B1==1] <- 2

  for (i in 1:nrow(B1)) {
    for (j in 1:ncol(B1)) {
      if (B1[i,j] == 1) B1[j,i] <- 2
      else if (B1[i,j] == 2) B1[j,i] <- 1
    }
  }

  # first two columns of pairs will contain indices of pairs
  # third column will contain whether the earlier is an effect of the later
  #    node (i.e. there is an edge, but no other possibly directed path in the
  #    PAG (i.e. directed path with tails or circles as endpoints) = 1,
  #    it is not a cause = 2,  we don't know = 3)
  # fourth column will contain what they are in the gen model (1=cause, ie.
  #    there is an edge or a path trough hidden variables, 2=not ancestor)
  pairs <- array(0,dim=c(nobs^2*(nrep-1),4))
  cnt <- 0

  # go through all the pairs with one variable in latest time layer and check
  # with true graph
  for (x in (nobs*(nrep-1)+1):(nobs*nrep)) { # later in time 
   
    for (y in 1:(nobs*(nrep-1))) { # earlier in time

      cnt <- cnt + 1
      pairs[cnt,1:2] <- c(y,x)

      if (pag[x,y] == 3) pairs[cnt,3] <- 3 # don't know
      if (pag[x,y] == 0) pairs[cnt,3] <- 2 # no direct effect
      if (pag[x,y] == 2) pairs[cnt,3] <- 2 # no direct effect   
      if (pag[x,y] == 1) { # potential cause
        # check if there is another possibly directed path from y to x, then we
        # don't know if direct cause
        if (isPossAncestorOf_indirect(y,x,pag,nrep)) pairs[cnt,3] <- 3 # don't know
        else pairs[cnt,3] <- 1 # direct effect
      }

      # check what it is in real graph
      if (isAncestorOf_cond(y,x,B1,nrep,1:(nrep*nobs)) ) {
        pairs[cnt,4] <- 1 # direct effect (possibly through hidden variables)
      }
      else pairs[cnt,4] <- 2 # not cause

    }
  }

  # give scores back:
  # 1st score: how often FCI could make a decision (is ancestor or not)
  # 2nd score: how often is this decision right (agrees with generating model)
  scores <- array(0,dim=c(1,2))

  cnt_possDec <- nrow(pairs)
  temp <- pairs[,3]!=3
  cnt_madeDec <- sum(temp)
  cnt_rightDec <- sum(pairs[temp,3]==pairs[temp,4])

#  print(pairs)

#  cat("possi decisions:", cnt_possDec, "\n")
#  cat("made  decisions:", cnt_madeDec, "\n")
#  cat("right decisions:", cnt_rightDec, "\n")

  scores[1,1] <- cnt_madeDec/cnt_possDec
  scores[1,2] <- cnt_rightDec/cnt_madeDec

  list(scores=scores,pair_res=pairs)

}



compareGenModelToPAG_ancestor <- function(B,pag,nrep,nobs) {

  # "ancestor"-score
  # B ... the connection matrix of the generating model over all variables
  # pag ... the PAG returned by (ts)(C)FCI over ALL observed variales

  # for each order pair (y,x) with y earlier in time than x, and x in last time
  # slice, check what decision we got from the PAG and compare it to the
  # generating model
  # 1 ... y ancestor of x (direct or indirect) = causal effect
  # 2 ... y not ancestor of x = no causal effect
  # 3 ... don't know (circles)

  # for comparison to generating model: put connection matrix nrep times in
  # one big matrix, similar to pag, with 1 and 2 as endpoints
  B1 <- rep_B_sort(t(B),nrep-1,nobs=nobs)
  B1[B1==1] <- 2

  for (i in 1:nrow(B1)) {
    for (j in 1:ncol(B1)) {
      if (B1[i,j] == 1) B1[j,i] <- 2
      else if (B1[i,j] == 2) B1[j,i] <- 1
    }
  }

  # first two columns of pairs will contain indices of pairs
  # third column will contain whether the earlier is a cause of the later
  #    node (i.e. there is an edge, but no other directed path in the PAG) = 1,
  #    it is not a cause = 2,  we don't know = 3
  # fourth column will contain what they are in the gen model (1=cause, ie.
  #    there is an edge or a path trough hidden variables, 2=not ancestor,
  #    0=didn't check)
  pairs <- array(0,dim=c(nobs^2*(nrep-1),4))
  cnt <- 0

  # go through all the pairs with one variable in latest time layer and check
  # with true graph
  for (x in (nobs*(nrep-1)+1):(nobs*nrep)) { # later in time 
   
    for (y in 1:(nobs*(nrep-1))) { # earlier in time

      cnt <- cnt + 1
      pairs[cnt,1:2] <- c(y,x)

      # get which decision the pag tells us, y earlier in time than x
      pathtype <- kindOfPath(x,y,pag,nrep,nobs)
      # if directed path, then ancestor (including directed edge)
      if (pathtype==1) pairs[cnt,3] <- 1
      # if partially directed path, then we don't know
      if (pathtype==3) pairs[cnt,3] <- 3
      # if all paths with at least one bi-directed edge, then non-cause
      # or if no path, then also non-cause
      if (pathtype==2) pairs[cnt,3] <- 2

      # check if y is ancestor of x in generating model
      if (isAncestorOf(y,x,B1,nrep)) pairs[cnt,4] <- 1
      else pairs[cnt,4] <- 2

    }
  }

  # give scores back:
  # 1st score: how often FCI could make a decision (is ancestor or not)
  # 2nd score: how often is this decision right (agrees with generating model)
  scores <- array(0,dim=c(1,2))

  cnt_possDec <- nrow(pairs)
  temp <- pairs[,3]!=3
  cnt_madeDec <- sum(temp)
  cnt_rightDec <- sum(pairs[temp,3]==pairs[temp,4])

#  print(pairs)

#  cat("possi decisions:", cnt_possDec, "\n")
#  cat("made  decisions:", cnt_madeDec, "\n")
#  cat("right decisions:", cnt_rightDec, "\n")

  scores[1,1] <- cnt_madeDec/cnt_possDec
  scores[1,2] <- cnt_rightDec/cnt_madeDec

  list(scores=scores,pair_res=pairs)

}



compareGenModelToPag_pairwise <- function(B, pag, nrep, nobs) {

  # pairwise score - compare each time-series variable (as in Granger)
  # B ... the connection matrix of the generating model over all variables
  # pag ... the PAG returned by (ts)(C)FCI over ALL observed variales

  # for each order pair (y,x) with x and y time-series variables (ie. collects
  # all lags of one variable: x = x(t-tau), ... , x(t)), check if y causes x
  # 1 ... y cause of x (with any lag) = causal effect
  # 2 ... y doesn't cause x (at any lag) = no causal effect
  # 3 ... don't know (circles)

  ntotal <- nobs*nrep

  # for comparison to generating model: put connection matrix nrep times in
  # one big matrix, similar to pag, with 1 and 2 as endpoints
  B1 <- rep_B_sort(t(B),nrep-1,nobs=nobs)
  B1[B1==1] <- 2

  for (i in 1:nrow(B1)) {
    for (j in 1:ncol(B1)) {
      if (B1[i,j] == 1) B1[j,i] <- 2
      else if (B1[i,j] == 2) B1[j,i] <- 1
    }
  }

  # real causal decisions from B
  # cause_B[y,x] = 1 iff y is a cause of x
  #              = 2 iff y s not a cause of x
  cause_B <- matrix(0,nobs,nobs)

  for (x in 1:nobs) {
    vx <- x + 0:(nrep-1)*nobs
    # only look at edges ending in the last variable, because of repeating str.
    xm <- max(vx)

    # check if y is a cause of x
    for (y in 1:nobs) {
      if (any(vx==y)) next
      vy <- y + 0:(nrep-1)*nobs
      # check what it is in real graph
      anc <- FALSE
      for (y0 in sort(vy,decreasing=TRUE)) {
        if (isAncestorOf_cond(y0,xm,B1,nrep,1:(nrep*nobs)) ) {
          cause_B[y,x] <- 1
          anc <- TRUE
          break
        }
      }
      if (!anc) { # not ancestor
        cause_B[y,x] <- 2
      }
    }
  }

  # causal decisions from PAG
  # cause_pag[y,x] = 1 iff y is a cause of x
  #                = 2 iff y is not a cause of x
  #                = 3 iff we don't know
  pag_save <- pag
  pag[upper.tri(pag)]<-0
  cause_pag <- matrix(0,nobs,nobs)

  for (x in 1:nobs) {
    vx <- x + 0:(nrep-1)*nobs
    # only look at edges ending in the last variable, because of repeating str.
    xm <- max(vx)

    # check if y is a cause of x
    for (y in 1:nobs) {
      if (any(vx==y)) next
      vy <- y + 0:(nrep-1)*nobs

      Pt <- pag[,vy][xm,] # contains all end-marks of edges pointing to xm

      if (any(Pt==1)) { # potential cause
        # check if there is another possibly directed path from y to x, then we
        # don't know if direct cause

        ind <- which(Pt==1)
        bool <- NULL
        for (i in ind) {
          bool <- c(bool,isPossAncestorOf_indirect(vy[i],xm,pag_save,nrep))
        }
        if (all(bool)) {
          cause_pag[y,x] <- 3 # don't know
        }
        else cause_pag[y,x] <- 1 # direct effect
      }
      else {
        if (any(Pt==3)) {
          cause_pag[y,x] <- 3
        }
        else {
          cause_pag[y,x] <- 2
        }
      }


    }

  }

  # give scores back:
  # 1st score: how often FCI could make a decision
  # 2nd score: how often is this decision right (agrees with generating model)
  scores <- array(0,dim=c(1,2))

  cnt_possDec <- nobs*(nobs-1)
  cnt_madeDec <- sum(cause_pag!=3) - nobs # diagonals are always 0
  cnt_rightDec <- sum(cause_pag==cause_B) - nobs # diagonals are always 0

#  cat("possi decisions:", cnt_possDec, "\n")
#  cat("made  decisions:", cnt_madeDec, "\n")
#  cat("right decisions:", cnt_rightDec, "\n")

  scores[1,1] <- cnt_madeDec/cnt_possDec
  scores[1,2] <- cnt_rightDec/cnt_madeDec

  list(cB = cause_B, cp = cause_pag, scores=scores)

}



compareGenModelToGranger_direct <- function(B,GG,nrep,nobs) {

  # "direct-cause" score
  # B ... the connection matrix of the generating model over all variables
  # GG ... the granger graph returned by Granger over ALL observed variales

  # for each order pair (y,x) with y earlier in time than x, and x in last time
  # slice, check what decision we got from the PAG and compare it to the
  # generating model
  # 1 ... y direct cause of x = causal effect
  # 2 ... y no direct cause of x = no causal effect
  # 3 ... don't know (circles)

  # for comparison to generating model: put connection matrix nrep times in
  # one big matrix, similar to pag, with 1 and 2 as endpoints
  B1 <- rep_B_sort(t(B),nrep-1,nobs=nobs)
  B1[B1==1] <- 2

  for (i in 1:nrow(B1)) {
    for (j in 1:ncol(B1)) {
      if (B1[i,j] == 1) B1[j,i] <- 2
      else if (B1[i,j] == 2) B1[j,i] <- 1
    }
  }

  # first two columns of pairs will contain indices of pairs
  # third column will contain whether the earlier is a direct cause of the later
  #    node = 1, it is not a cause = 2, 
  # fourth column will contain what they are in the gen model (1=cause, ie.
  #    there is an edge or a path trough hidden variables, 2=not ancestor)
  pairs <- array(0,dim=c(nobs^2*(nrep-1),4))
  cnt <- 0

  # go through all the pairs with one variable in latest time layer and check
  # with true graph
  for (x in (nobs*(nrep-1)+1):(nobs*nrep)) { # later in time 
   
    for (y in 1:(nobs*(nrep-1))) { # earlier in time

      cnt <- cnt + 1
      pairs[cnt,1:2] <- c(y,x)

      if (GG[x,y] == 1) pairs[cnt,3] <- 1 # direct effect
      else pairs[cnt,3] <- 2 # no direct effect   

      # check what it is in real graph
      if (isAncestorOf_cond(y,x,B1,nrep,1:(nrep*nobs)) ) {
        pairs[cnt,4] <- 1 # direct effect (possibly through hidden variables)
      }
      else pairs[cnt,4] <- 2 # not cause

    }
  }

  # give scores back:
  # 1st score: how often FCI could make a decision (is ancestor or not)
  # 2nd score: how often is this decision right (agrees with generating model)
  scores <- array(0,dim=c(1,2))

  cnt_possDec <- nrow(pairs)
  temp <- pairs[,3]!=3
  cnt_madeDec <- sum(temp)
  cnt_rightDec <- sum(pairs[temp,3]==pairs[temp,4])

#  print(pairs)

#  cat("possi decisions:", cnt_possDec, "\n")
#  cat("made  decisions:", cnt_madeDec, "\n")
#  cat("right decisions:", cnt_rightDec, "\n")

  scores[1,1] <- cnt_madeDec/cnt_possDec
  scores[1,2] <- cnt_rightDec/cnt_madeDec

  list(scores=scores,pair_res=pairs)

}



compareGenModelToGranger_ancestor <- function(B,GG,nrep,nobs) {

  # "ancestor"-score
  # B ... the connection matrix of the generating model over all variables
  # GG ... the granger graph returned by Granger over ALL observed variales

  # for each order pair (y,x) with y earlier in time than x, and x in last time
  # slice, check what decision we got from the PAG and compare it to the
  # generating model
  # 1 ... y ancestor of x (direct or indirect) = causal effect
  # 2 ... y not ancestor of x = no causal effect
  # 3 ... don't know (circles)

  # for comparison to generating model: put connection matrix nrep times in
  # one big matrix, similar to pag, with 1 and 2 as endpoints
  B1 <- rep_B_sort(t(B),nrep-1,nobs=nobs)
  B1[B1==1] <- 2

  for (i in 1:nrow(B1)) {
    for (j in 1:ncol(B1)) {
      if (B1[i,j] == 1) B1[j,i] <- 2
      else if (B1[i,j] == 2) B1[j,i] <- 1
    }
  }

  # first two columns of pairs will contain indices of pairs
  # third column will contain whether the earlier is an ancestor of the later
  #    node = 1 (ie. directed path), it is not an ancestor = 2
  # fourth column will contain what they are in the gen model (1=cause, ie.
  #    there is an edge or a path trough hidden variables, 2=not ancestor)
  pairs <- array(0,dim=c(nobs^2*(nrep-1),4))
  cnt <- 0

  # go through all the pairs with one variable in latest time layer and check
  # with true graph
  for (x in (nobs*(nrep-1)+1):(nobs*nrep)) { # later in time 
   
    for (y in 1:(nobs*(nrep-1))) { # earlier in time

      cnt <- cnt + 1
      pairs[cnt,1:2] <- c(y,x)

      # get which decision the pag tells us, y earlier in time than x
      pathtype <- kindOfPath(x,y,GG,nrep,nobs)
      # if directed path, then ancestor (including directed edge)
      if (pathtype==1) pairs[cnt,3] <- 1
      # if all paths with at least one bi-directed edge, then non-cause
      # or if no path, then also non-cause
      if (pathtype==2) pairs[cnt,3] <- 2

      # check if y is ancestor of x in generating model
      if (isAncestorOf(y,x,B1,nrep)) pairs[cnt,4] <- 1
      else pairs[cnt,4] <- 2

    }
  }

  # give scores back:
  # 1st score: how often FCI could make a decision (is ancestor or not)
  # 2nd score: how often is this decision right (agrees with generating model)
  scores <- array(0,dim=c(1,2))

  cnt_possDec <- nrow(pairs)
  temp <- pairs[,3]!=3
  cnt_madeDec <- sum(temp)
  cnt_rightDec <- sum(pairs[temp,3]==pairs[temp,4])

#  print(pairs)

#  cat("possi decisions:", cnt_possDec, "\n")
#  cat("made  decisions:", cnt_madeDec, "\n")
#  cat("right decisions:", cnt_rightDec, "\n")

  scores[1,1] <- cnt_madeDec/cnt_possDec
  scores[1,2] <- cnt_rightDec/cnt_madeDec

  list(scores=scores,pair_res=pairs)

}



compareGenModelToGranger_pairwise <- function(B,GG,nrep,nobs) {

  # pairwise score - compare each time-series variable (as in Granger)
  # B ... the connection matrix of the generating model over all variables
  # GG ... the granger graph returned by Granger over ALL observed variales

  # for each order pair (y,x) with x and y time-series variables (ie. collects
  # all lags of one variable: x = x(t-tau), ... , x(t)), check if y causes x
  # 1 ... y cause of x (with any lag) = causal effect
  # 2 ... y doesn't cause x (at any lag) = no causal effect
  # 3 ... don't know (circles)


  ntotal <- nobs*nrep

  # for comparison to generating model: put connection matrix nrep times in
  # one big matrix, similar to pag, with 1 and 2 as endpoints
  B1 <- rep_B_sort(t(B),nrep-1,nobs=nobs)
  B1[B1==1] <- 2

  for (i in 1:nrow(B1)) {
    for (j in 1:ncol(B1)) {
      if (B1[i,j] == 1) B1[j,i] <- 2
      else if (B1[i,j] == 2) B1[j,i] <- 1
    }
  }

  # real causal decisions from B
  # cause_B[y,x] = 1 iff y is a cause of x
  #              = 2 iff y s not a cause of x
  cause_B <- matrix(0,nobs,nobs)

  for (x in 1:nobs) {
    vx <- x + 0:(nrep-1)*nobs
    # only look at edges ending in the last variable, because of repeating str.
    xm <- max(vx)

    # check if y is a cause of x
    for (y in 1:nobs) {
      if (any(vx==y)) next
      vy <- y + 0:(nrep-1)*nobs
      # check what it is in real graph
      anc <- FALSE
      for (y0 in sort(vy,decreasing=TRUE)) {
        if (isAncestorOf_cond(y0,xm,B1,nrep,1:(nrep*nobs)) ) {
          cause_B[y,x] <- 1
          anc <- TRUE
          break
        }
      }
      if (!anc) {
        cause_B[y,x] <- 2
      }
    }
  }

  # causal decisions from GG
  # cause_GG[y,x] = 1 iff y is a cause of x
  #                = 2 iff y is not a cause of x
  GG_save <- GG
  GG[upper.tri(GG)]<-0
  cause_GG <- matrix(0,nobs,nobs)

  for (x in 1:nobs) {
    vx <- x + 0:(nrep-1)*nobs
    # only look at edges ending in the last variable, because of repeating str.
    xm <- max(vx)

    # check if y is a cause of x
    for (y in 1:nobs) {
      if (any(vx==y)) next
      vy <- y + 0:(nrep-1)*nobs

      Gt <- GG[,vy][xm,] # contains all end-marks of edges pointing to xm

      if (any(Gt==1)) { 
        cause_GG[y,x] <- 1
      }
      else {
        cause_GG[y,x] <- 2
      }

    }

  }

  # give scores back:
  # 1st score: how often could make a decision - is one
  # 2nd score: how often is this decision right (agrees with generating model)
  scores <- array(0,dim=c(1,2))

  cnt_possDec <- nobs*(nobs-1)
  cnt_madeDec <- sum(cause_GG!=3) - nobs # diagonals are always 0
  cnt_rightDec <- sum(cause_GG==cause_B) - nobs # diagonals are always 0

#  cat("possi decisions:", cnt_possDec, "\n")
#  cat("made  decisions:", cnt_madeDec, "\n")
#  cat("right decisions:", cnt_rightDec, "\n")

  scores[1,1] <- cnt_madeDec/cnt_possDec
  scores[1,2] <- cnt_rightDec/cnt_madeDec

  list(cB = cause_B, cG = cause_GG, scores=scores)

}
