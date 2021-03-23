# This file contains some useful functions for time series graphs (at the 
# moment only 1st order)
# 
# input:
# B ... matrix containing the effects from time t-1 to time t
# for more input details see the specific functions
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


parents <- function(B,Z) {
  # gives back the parents of the set Z in B
  # Z is a vector
  # if no parents, return NULL
  len <- length(Z)
  pa.Z <- NULL
  for (i in 1:len) {
    pa.Z <- c(pa.Z,which(B[Z[i],]==1))
  }
  sort(unique(pa.Z))
}


rep_B <- function(B, nrep) {
  # gives back the DAG-matirx for a time series repeated nrep times
  # nrep: how often _edges_ are repeated, ie nodes are repeated nrep+1 times

  nvar <- ncol(B)

  Brep <- kronecker(diag(nrep),B)
  Brep <- cbind(array(0,dim=c(nrep*nvar,nvar)), Brep)
  Brep <- rbind(Brep, array(0,dim=c(nvar,nrep*nvar+nvar)))
  Brep
}


rep_B_sort <- function(B, nrep, nobs=2) {
  # similar to rep_B but puts first all observed and then all hidden variables
  # nrep: how often _edges_ are repeated, ie nodes are repeated nrep+1 times
  # NOTE: this only works for models with 1 lag!

  nvar <- ncol(B)
  nhid <- nvar-nobs

  Brep <- array(0,dim=c(nvar*nrep,nvar*nrep))
  Brep <- kronecker(diag(nrep),B)

  if (nhid == 0) {
    Brep <- cbind(array(0,dim=c(nrep*nvar,nobs)),Brep)
    Brep <- rbind(Brep,array(0,dim=c(nobs,ncol(Brep))))
    return(Brep)
  }

  # need to row- and column-permute Brep, such that observed variables come
  # before hidden variables and include some empty rows and columns
  vec_obs <- 1:nobs
  vec <- NULL
  for (i in 1:nrep) vec <- c(vec,vec_obs+nvar*(i-1))
  vec_diff <- setdiff(1:(nvar*nrep),vec)
  vec <- c(vec, vec_diff) # new order of variables

  perm <- array(0,dim=c(nvar*nrep,nvar*nrep))
  for (i in 1:(nvar*nrep)) perm[vec[i],i] <- 1
  Brep <- t(perm)%*%Brep%*%perm
  Brep <- cbind(array(0,dim=c(nrep*nvar,nobs)), # no edges ouo of 1st obs. var
                Brep[,1:(nobs*nrep)],
                array(0,dim=c(nrep*nvar,nhid)), # no edges out of 1st hid. var
                Brep[,(nobs*nrep+1):(nrep*nvar)])
  Brep <- rbind(Brep[1:(nobs*nrep),],
                array(0,dim=c(nobs,ncol(Brep))), # no edges into last obs. var
                Brep[(nobs*nrep+1):(nvar*nrep),],
                array(0,dim=c(nhid,ncol(Brep)))) # no edges into last hid. var
  class(Brep) <- 'hDAG'
  Brep
}


common_ancestor <- function(B,x,y) {

  # returns true if nodes x and y have a common ancestor in B, looking
  # n(n-1)/2 steps (edges) back (n=nvar)
  # x and y in same time slice, 1-lag model. Only check if when going back one
  # step, if the ancestorsets overlap

  nvar <- nrow(B)

  anc1 <- parents(B,x)
  anc2 <- parents(B,y)

  for (k in 1:(nvar*(nvar-1)/2)) { # how often _edges_ are repeated

    if (length(intersect(anc1,anc2))!=0) { # if ancestors intersect
      return(TRUE) # found a common ancestor
    }

    anc1 <- parents(B,anc1) # only need ancestors in time layer k
    anc2 <- parents(B,anc2)

  } # end k

  return(FALSE)

}


get_closure <- function(B) {

  nvar <- ncol(B)

  closure <- array(0,dim=c(nvar,0))
  for (m in 1:(nvar-1)) {
    for (n in (m+1):nvar) {
      if (common_ancestor(B,m,n)) {
        # create new variable connecting node m and n
        temp <- array(0,dim=c(nvar,1))
        temp[c(m,n)] <- 1
        closure <- cbind(closure,temp)
      }
    }
  }

  closure

}


add_closure <- function(Brep,closure,nvar,nobs,twindow) {

  # Brep: a ts-DAG-matrix, for example from rep_B_sort
  # closure: a matrix with nvar rows and at least 1 column telling which nodes
  #   have a common ancestor, f.ex. from get_closure
  # nvar: total number of variables
  # nobs: observed number of variables
  # twindow: how often _nodes_ are repeated in Brep

  nvar_all <- nvar*twindow
  nobs_all <- nobs*twindow
  nhid <- nvar-nobs

  closure_aug <- array(0,dim=c(nvar_all,ncol(closure)))
  closure_aug[(nobs_all-nobs+1):nobs_all,] <- closure[1:nobs,]
  if (nvar_all!=nobs_all) {
    closure_aug[(nvar_all-nhid+1):nvar_all,] <- closure[(nobs+1):nvar,]
  }

  Brep <- cbind(Brep,closure_aug)
  Brep <- rbind(Brep,array(0,dim=c(ncol(closure),ncol(Brep))))

  Brep

}


add_closure_new <- function(Brep,closure,nvar,nobs,twindow) {

  # Brep: a ts-DAG-matrix, for example from rep_B_sort
  # closure: a matrix with nvar rows and at least 1 column telling which nodes
  #   have a common ancestor, f.ex. from get_closure
  # nvar: total number of variables
  # nobs: observed number of variables
  # twindow: how often _nodes_ are repeated in Brep

  nvar_all <- nvar*twindow
  nobs_all <- nobs*twindow
  nhid <- nvar-nobs

  closure_aug <- array(0,dim=c(nvar_all,ncol(closure)))
  closure_aug[1:nobs,] <- closure[1:nobs,]
  if (nvar_all!=nobs_all) {
    closure_aug[(nobs_all+1):(nobs_all+nhid),] <- closure[(nobs+1):nvar,]
  }

  Brep <- cbind(Brep,closure_aug)
  Brep <- rbind(Brep,array(0,dim=c(ncol(closure),ncol(Brep))))

  Brep

}



haveConfounder <- function(x,y,B,vec,nrep) {

  # check if x and y (being in arbitrary time slice) have a common ancestor in 
  # B given a vector vec of conditioning nodes (if vec lies on path to ancestor 
  # this path is blocked)

  anc_x <- setdiff(parents(B,x),vec)
  anc_y <- setdiff(parents(B,y),vec)

  for (i in 1:nrep) {
    anc_x <- union(anc_x,setdiff(parents(B,anc_x),vec))
    anc_y <- union(anc_y,setdiff(parents(B,anc_y),vec))
  }
  if(length(intersect(anc_x, anc_y)) > 0) return(TRUE)
  else return(FALSE)

}


ancestors <- function(B, Z, anc=Z) {
  # returns the ancestors of a set Z (NB! Z is itself an ancestor too.)
  # if no ancestors, then returns Z

  pa.Z <- NULL
  for (i in 1:length(Z)) {
    pa.Zi <- which(B[Z[i],]==1)
    pa.Z <- c(pa.Z, pa.Zi)
  }

  anc <- union(anc,pa.Z)

  if (length(pa.Z)!=0) ancestors(B, pa.Z, anc)
  else return(anc)

} # end fct


isAncestorOf_cond <- function(x,y,B,nrep,vec=NULL) {

  # check if x is an ancestor of y in B given a vector vec of cond. nodes (ie.
  # can't go through this nodes when looking for ancestors)

  anc_y <- parents(B,y)

  for (i in 1:nrep) {
    if (length(anc_y)==0) return(FALSE)
    if (is.element(x,anc_y)) return(TRUE)
    else {
      anc_y <- setdiff(anc_y,vec)
      anc_y <- parents(B,anc_y)
    }
  }

  FALSE

}


isAncestorOf <- function(x,y,B,nrep) {

  # check if x is an ancestor of y in B

  anc_y <- which(B[y,]==1)

  for (i in 1:nrep) {
    if (length(anc_y)==0) return(FALSE)
    if (is.element(x,anc_y)) return(TRUE)
    else {
      temp <- NULL
      for (j in 1:length(anc_y)) {
        temp <- union(temp,which(B[anc_y[j],]==1))
      }
      anc_y <- temp
    }
  }

  FALSE # is not ancestor

}


isAncestorOf_indirect <- function(x,y,B,nrep) {

  # check if x in an indirect ancestor of y (possible also direct, but want to
  # find path of length at least 2

  anc_y <- which(B[y,]==1)
  # if x is now in anc_y then there is a direct edge from x to y, don't want to
  # consider this, take x out of anc_y
  anc_y <- setdiff(anc_y,x)

  for (i in 1:nrep) {
    if (length(anc_y)==0) return(FALSE)
    if (is.element(x,anc_y)) return(TRUE)
    else {
      temp <- NULL
      for (j in 1:length(anc_y)) {
        temp <- union(temp,which(B[anc_y[j],]==1))
      }
      anc_y <- temp
    }
  }

  FALSE # is not ancestor
}


isPossAncestorOf_indirect <- function(x,y,B,nrep) {

  # check if x is a possible indirect ancestor of y (possible also direct, but
  # want to find path of length at least 2), that means there is a partially
  # directed path from x to y of length at least 2

  anc_y <- which(B[y,]==1 | B[y,]==3)
  # if x is now in anc_y then there is a direct edge from x to y, don't want to
  # consider this, take x out of anc_y
  anc_y <- setdiff(anc_y,x)

  for (i in 1:nrep) {
    if (length(anc_y)==0) return(FALSE)
    if (is.element(x,anc_y)) return(TRUE)
    else {
      temp <- NULL
      for (j in 1:length(anc_y)) {
        temp <- union(temp,which(B[anc_y[j],]==1 | B[anc_y[j],]==3))
      }
      anc_y <- temp
    }
  }

  FALSE # is not ancestor
}


getTier <- function(x,nrep,nobs) {

  # gives back the index of the tier of a variable
  # x must be a number between 1 and nrep*nobs, since these are all obs. var.

  tier <- max(which(0:(nrep-1)*nobs+1 <= x))

  tier

}


kindOfPath <- function(x,y,pag,nrep,nobs) {
  # what kind of path is from x (later in time) back to y (earlier in time)
  # 1 ... directed path: y-->a-->...-->b-->x
  #       y causes x
  # 2 ... "bidirected" path: yo->a<->...-->b<->x (at least one bidirected edge
  #       per path) or no path at all
  #       y is not a cause of x
  # 3 ... parially directed path: yo->a-->...o->b-->x (only o-> and -->)
  #       we don't know

  # Check first for a directed path, because if this is present, then we always
  # can say that y causes x. Then check for a partially directed paths, because
  # if this is present, we can never know if it would be a direct cause or not.
  # If there is neither a directed path nor a parthially directed path, we know
  # that y is not a cause of x. (? cover really all the cases?)

  ntier_x <- getTier(x,nrep,nobs)

  # check if there is a directed path y-->a-->...-->b-->x

  anc_x <- which(pag[x,]==1) # set of nodes with z-->x

  # check if y is an ancestor of x (could have used function...)
  for (i in 1:nrep) {
    if (length(anc_x) == 0) break;
    if (is.element(y,anc_x)) return(1) # i.e. there is directed path
    else { # find ancestors of ancestors
      temp <- NULL
      for (j in anc_x) {
        temp <- union(temp,which(pag[j,]==1))
      }
      anc_x <- temp
    }
  }

  # if y is not an ancestor of x, check if there is a partially directed path 
  anc_x <- which(pag[x,]==1)
  cir_x <- which(pag[x,]==3) # nodes with zo->x

  anc_cir_x <- union(anc_x,cir_x)

  for (i in 1:nrep) {
    if (length(anc_cir_x) == 0) break;
    if (is.element(y,anc_cir_x)) return(3) # i.e. there is partially dir. path
    else { # go further back
      temp <- NULL
      temp1 <- NULL
      for (j in anc_cir_x) {
        temp <- union(temp,which(pag[j,]==1))
        temp1 <- union(temp1,which(pag[j,]==3))
      }
      anc_cir_x <- union(temp,temp1)
    }
  }

  # if y is not an ancestor of x, and we are not uncertain, then y is not a
  # cause of x  

  return(2)

}


findMarkovBlanket <- function(x,B,nobs,nvar,nrep) {

  # finds the Markov Blanket of node x in tsDAG B
  # is specialized to case where x is in latest timer layer

  nhid <- nvar-nobs
  nhid_all <- nhid*nrep
  nobs_all <- nobs*nrep

  possMB <- 1:(nobs*(nrep-1))
  possMB <- union(possMB,(nobs*nrep+1):(nobs*nrep+nhid_all-nhid))

  notMB <- setdiff(1:(nvar*nrep),possMB)

  MB <- NULL
  inMBobs <- NULL
  inMBhid <- NULL

  pa.x <- which(B[x,]==1)
  ch.indic <- rep(0,length(pa.x)) # has a one if corresponding node is a child
                                  # ie arrow is into it

  if (length(pa.x)==0) return(NULL)
  flag <- TRUE

  while (flag) {

    y <- pa.x[1]
    ind <- ch.indic[1]
    pa.x <- pa.x[index(2,length(pa.x))]
    ch.indic <- ch.indic[index(2,length(ch.indic))]

    if (is.element(y,possMB) | ind==1) {
      if (dconnected(B, x, y, Z=inMBobs)) {
        MB <- c(MB,y)

        if (y <= nobs*(nrep-1)) {
          inMBobs <- c(inMBobs,y)
          if (ind==1) {
            pa.y <- setdiff(which(B[y,]==1),union(pa.x,MB))
            pa.x <- c(pa.x,pa.y)
            ch.indic <- c(ch.indic,rep(0,length(pa.y)))
          }
        }
        else {
          inMBhid <- c(inMBhid,y)
          if (ind!=1) {
            pa.y <- setdiff(which(B[y,]==1),union(pa.x,MB))
          }
          ch.y <- setdiff(which(B[,y]==1),union(union(pa.x,notMB),inMBhid))
          pa.x <- c(pa.x,pa.y,ch.y)
          ch.indic <- c(ch.indic,rep(0,length(pa.y)),rep(1,length(ch.y)))
        }

        possMB <- setdiff(possMB,y) # don't have to go there again
      }
    }

    if (length(pa.x)==0) flag <- FALSE

  }

  return(sort(unique(MB)))

}


index<-function(from,to) {
  if ( from > to ) {
    R<-c()
  } else {
    R<-from:to
  }
  R
}
