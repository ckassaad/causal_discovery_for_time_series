# This function returns TRUE if x and y are dconnected given Z
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



dconnected <- function(B, x, y, Z=NULL, income.x=FALSE, path=NULL,
                       ancestorZ=NULL) {

  if (is.null(ancestorZ)) ancestorZ <- ancestors(B,Z)

  # reached x from y through a d-connecting path
  if (x == y) {
    return(TRUE)
  }

  # visited x already and couldn't find d-connecting path
  if (x %in% path) { 
    return(FALSE)
  }

  path <- c(path,x)

  pa.x <- which(B[x,]==1)
  ch.x <- which(B[,x]==1)

  edges.x <- array(0,dim=c(length(pa.x)+length(ch.x),2))
  edges.x[,1] <- c(pa.x,ch.x)
  edges.x[,2] <- c( rep(1,length(pa.x)), rep(0,length(ch.x)))

  for (i in index(1,nrow(edges.x))) {
    if (edges.x[i,2] == 1) income.from.neigh <- TRUE
      # where TRUE means that edge is pointing towards x
    else income.from.neigh <- FALSE

    # check if there is a collider on the path: comefrom -> x <- neighbour
    isCollider <- income.x & income.from.neigh
    # if there is collider, check if path is d-connected
    passAsCollider <- isCollider & (x %in% ancestorZ)
    # if there is no collider, check if path is d-connected
    passAsNonCollider <- !isCollider & (!x %in% Z)

    if (passAsCollider | passAsNonCollider) { # ie. dconnected so far
      new.x <- edges.x[i,1] # second node in edge besides x
      # get endpoint of edge: f. ex. if x -> new.x, than income.new.x = TRUE
      if (edges.x[i,2] == 1) income.new.x <- FALSE 
      else income.new.x <- TRUE

      if (dconnected(B, new.x, y, Z, income.new.x, path, ancestorZ)) {
        return(TRUE)
      }

    } # end if (passAsCollider | passAsNonCollider)

  } #end "for i"

  path <- path[index(1,length(path)-1)]
  return(FALSE)

}
