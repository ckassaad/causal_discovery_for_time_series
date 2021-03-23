# This function plots a time series graph
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


plot_timeseries <- function(lis, nrepeat=0, nodelabels=NULL,
    ngvec=rep(FALSE,ncol(lis[[1]])), filename='tempfile',
    hidden=FALSE, nobs=NULL) {

  # lis ... list with:
  #          first element = matrix of instantaneous effects
  #          second element = matrix with 1 time lag
  #          third element = matrix with 2 time lags
  #          and so on

  # nobs ... number of observed variable - only temporary input, should be
  #          automatically calculated from lis, only needed when hidden vars
  #          are also included

  # nrepeat ... how often each lag shall be additionally plotted (repeated)
  #             (plot the same time series once plus "nrepeat" times after each
  #             other

  # nodelabels ... list of strings given the names of variables
  #                or it NULL, then using nodelabels "1", "2", ...
  # example: lis[[1]]=B0, lis[[2]]=B1, B0,B1 4by4 matrices
  #          nodelabels[[1]]="empl",nodelabels[[2]]="sales",
  #          nodelabels[[3]]="rnd",nodelabels[[4]]="opinc"

  # ngvec ... boolean, indicate which variables are observed (FALSE) and which
  #           are hidden (TRUE)

  # filename ... the name where the plot gets saved (filename.dot, filename.ps)

  # hidden ... boolean, indicates whether hidden variables are included (TRUE)
  #            or not (FALSE)

  nrepeat <- 1+nrepeat # plot it once as usuall and nrepeat times more

  len <- length(lis)
  for (i in 1:len) {
     if ( all(dim(lis[[1]]) != dim(lis[[i]])) ) {
       return(cat("Dimensions of coeffient matrices don\'t agree!\n"))
     }
#     if ( !is.null(nodelabels) & length(nodelabels) != dim(lis[[1]])[1] ) {
#       return(cat("Not enough or too many nodelabels!\n"))
#     }
  }

  if (!hidden & is.null(nobs)) nobs <- dim(lis[[1]])[1]

  nlags <- len-1 # number of lags
  nhid <- dim(lis[[1]])[1] - nobs
  ntotal <- nobs + nhid

  if ( !is.null(nodelabels) ) {
    if (length(nodelabels) == dim(lis[[1]])[1]) {
      nodelabels <- rep(nodelabels,nrepeat+nlags)
    }
    else {
      if (length(nodelabels) != dim(lis[[1]])[1]*(nrepeat+nlags)) {
        return(cat("Not enough or too many nodelabels!\n"))
      }
    }
  }

  if (hidden) {
    # cut out hidden variables which aren't connected to any observed variable
    # in any time lag (don't plot those)
    indcut <- NULL
    for (i in (nobs+1):ntotal) {
      temp <- TRUE
      for (j in 0:nlags) {
        if( !all(lis[[j+1]][,i]==0) | !all(lis[[j+1]][i,]==0) ) {
          temp <- FALSE
          break
        }
      }
      if (temp) indcut <- c(indcut,i)
    }
    for (j in 0:nlags) {
      ntake <- setdiff(1:ntotal,indcut)
      lis[[j+1]] <- lis[[j+1]][ntake,ntake]
    }

    # indicate how many observed and hidden variables for node-shapes
    ngvec <- c(rep(F,nobs),rep(T,dim(lis[[1]])[1]-nobs))
  }

  # generate the dot file
  dotfilename <- paste(sep='',filename,'.dot')
  psfilename <- paste(sep='',filename,'.ps')

  # to distinguish between observed and latent variables
  nodeshapevec <- ngvec

  # if TRUE arcs are labeled with the strength of connection (doesn't look that
  # nice yet, often numbers overlap...)
  usearclabels <- FALSE #TRUE


  # writes nodes and adjacencies in file dotfilename(tempfile.dot)
  timeseries2dot(lis, len, nrepeat, nodelabels, nodeshapevec,
    usearclabels, dotfilename)


  # make a ps file using graphviz
  useprog <- 'neato' #-s -n2
  progcall <- paste(sep='',useprog,' -Tps -s72 ',dotfilename,' -o ',psfilename)
  system(progcall)

  # show the ps file using gv
  gvcall <- paste(sep='','gv ',psfilename)
  system(gvcall,wait=FALSE)

}



plot_timeseries_closure <- function(lis, closure, nrepeat=0, nodelabels=NULL,
    ngvec=rep(1,ncol(lis[[1]])), filename='tempfile', hidden=FALSE,
    nobs=NULL) {

  # inputs: see plot_timeseries
  # additionally:
  #   closure ... matrix for which variables have a common ancestor, to
  #               be replaced with a dummy (closure) variable, get f.ex.
  #               with function get_closure inf ts_functions.R
  #   ngvec ... 3 values: 1 -> observed variable, round node
  #                       2 -> hidden variable, squared node
  #                       3 -> closure variable, diamond node

  nrepeat <- 1+nrepeat # plot it once as usuall and nrepeat times more

  len <- length(lis)
  for (i in 1:len) {
     if ( all(dim(lis[[1]]) != dim(lis[[i]])) ) {
       return(cat("Dimensions of coeffient matrices don\'t agree!\n"))
     }
     if ( !is.null(nodelabels) & length(nodelabels) != dim(lis[[1]])[1] ) {
       return(cat("Not enough or too many nodelabels!\n"))
     }
  }

  if (!hidden & is.null(nobs)) nobs <- dim(lis[[1]])[1]

  # if (hidden): need nobs given as input parameter
  nlags <- len-1 # number of lags
  nvar <- dim(lis[[1]])[1]
  nhid <- nvar - nobs
  nclo <- ncol(closure) # number of closure-variables

  # indicate how many observed and hidden variables for node-shapes
  ngvec <- c(rep(1,nobs),rep(2,nhid),rep(3,nclo))

  if (hidden) {
    nodelabels <- NULL
    for (i in 1:(nlags+nrepeat)) {
      temp_hv <- NULL
      for (j in nhid:1) {
        temp_hv <- c(temp_hv,paste(sep='','hv',(i-1)*nhid+j))
      }
#      nodelabels <- c(nodelabels, c(temp_hv, nobs:1+(i-1)*nobs))
      nodelabels <- c(c(temp_hv, nobs:1+(i-1)*nobs),nodelabels)
    }
  }

  # generate the dot file
  dotfilename <- paste(sep='',filename,'.dot')
  psfilename <- paste(sep='',filename,'.ps')

  # to distinguish between observed, latent and closure variables
  nodeshapevec <- ngvec

  # if TRUE arcs are labeled with the strength of connection (doesn't look that
  # nice yet, often numbers overlap...)
  usearclabels <- FALSE #TRUE


  # writes nodes and adjacencies in file dotfilename(tempfile.dot)
  timeseries_closure2dot(lis, closure, len, nrepeat, nodelabels, nodeshapevec,
    usearclabels, dotfilename)


  # make a ps file using graphviz
  useprog <- 'neato' #-s -n2
  progcall <- paste(sep='',useprog,' -Tps -s72 ',dotfilename,' -o ',psfilename)
  system(progcall)

  # show the ps file using gv
  gvcall <- paste(sep='','gv ',psfilename)
  system(gvcall,wait=FALSE)

}



timeseries2dot <- function(lis, len, nrepeat, nodelabels, nodeshapevec, 
    usearclabels, filename) {

  # number of nodes in time series (per time step)
  nnodes <- dim(lis[[1]])[1]

  # number of time lags in the time series
  nlags <- len - 1

  # if no nodelabels are given, create them as numbers
  if (is.null(nodelabels)){
    nodelabels <- list()
    for (i in 1:nnodes) {
      nodelabels[[i]] <- sprintf("%d",i)
    }
    nodelabels <- rep(nodelabels,nrepeat+nlags)
  }

  # node shapes: boolean vector indicating FALSE=round for observed variables
  # or TRUE=square for latent variables
  nodeshapes <- list()
  for (i in 1:nnodes) {
    if (nodeshapevec[i]) {
     nodeshapes[[i]] <- 'box'
    }
    else {
      nodeshapes[[i]] <- 'ellipse'
    }
  }

  # usearclabels is boolean, if TRUE we put numerical values on arcs
  if (usearclabels) {
    arclabels <- list()
    for (l in 1:len) {
      for (i in 1:(nnodes*nnodes)) {
        if (l==1) arclabels[[i]] <- sprintf("%.3f",lis[[l]][i])
        else arclabels[[i]][l] <- sprintf("%.3f",lis[[l]][i])
      }
    }
  }

  # construct a format string for nodes
  nodeformat <- ' %d [pos = "%f,%f", label = "%s", shape = %s, style = "%s",
    pin = true];'

  # construct a format string for edges
  if (usearclabels) {
    attributes <- 'style = %s, color = %s, label = "%s", dir = %s'
  }
  else {
    attributes <- 'style = %s, color = %s, dir = %s'
  }
  edgeformat = paste(sep='','  %d -> %d [', attributes, '];')

  # write the beginning of the file
  write('digraph G {', file=filename, append=FALSE, sep="")
  write('  center = 1;', file=filename, append=TRUE, sep="")
  write('  size = \"10, 10\";', file=filename, append=TRUE, sep="")
  write('  splines = true;', file=filename, append=TRUE, sep="")
  write('  sep = 1;', file=filename, append=TRUE, sep="")

  # write the nodes
  for (l in 1:(nrepeat+nlags)) {
    for (node in 1:nnodes) {
      thisnode <- sprintf(nodeformat, node+nnodes*(l-1), 100+200*(l-1),
                  100*(nnodes-node), nodelabels[[node+(l-1)*(nnodes)]],
                  nodeshapes[[node]], 'solid')
      write(thisnode, file=filename, append=TRUE, sep="")
    }
  }


  # write the edges
  threshold <- 0 # plot all edges 
             # 0.1 # only edges with connection absolute stronger than 0.1

  for (l in 1:len) {
    for (r in 1:(nrepeat+nlags-l+1)) {
      for (node1 in 1:nnodes) {
        for (node2 in 1:nnodes) {
          if (abs(lis[[l]][node1,node2]) > threshold) {
            #style <- 'solid'
            if(lis[[l]][node1,node2] > 0) {
              color <- 'black' #'black' #'green'
              style <- 'solid'
            }
            else {
              color <- 'black' #'black' #'red'
              style <- 'dashed'
            }
            if (usearclabels) {
              thisarc <- sprintf(edgeformat, node2+(r-1)*nnodes,
                     node1+nnodes*(l-1)+(r-1)*nnodes, style, color,
                     arclabels[[(node2-1)*nnodes+node1]][l], 'forward')
            }
            else {
              thisarc <- sprintf(edgeformat, node2+(r-1)*nnodes,
                     node1+nnodes*(l-1)+(r-1)*nnodes, style, color, 'forward')
             }
             write(thisarc, file=filename, append=TRUE, sep="")
          } # end if (abs(...)>threshold)
        } # end node2
      } # end node1
    } # end r
  } # end l

  write('}', file=filename, append=TRUE, sep="")

}



timeseries_closure2dot <- function(lis, closure, len, nrepeat, nodelabels,
    nodeshapevec, usearclabels, filename) {

  # number of nodes in time series (per time step)
  nnodes <- dim(lis[[1]])[1]

  # number of time lags in the time series
  nlags <- len - 1
  nclo <- ncol(closure)

  # if no nodelabels are given, create them as numbers
  if (is.null(nodelabels)){
    nodelabels <- list()
    for (i in 1:nnodes) {
      nodelabels[[i]] <- sprintf("%d",i)
    }
    nodelabels <- rep(nodelabels,nrepeat+nlags)
  }

  # node shapes: 1 = observed var, 2 = hidden var, 3 = closure var
  nodeshapes <- list()
  for (i in 1:(nnodes+nclo)) {
    if (nodeshapevec[i] == 1) {
      nodeshapes[[i]] <- 'ellipse'
    }
    if (nodeshapevec[i] == 2) {
      nodeshapes[[i]] <- 'box'
    }
    if (nodeshapevec[i] == 3) {
      nodeshapes[[i]] <- 'diamond'
    }
  }

  # usearclabels is boolean, if TRUE we put numerical values on arcs
  if (usearclabels) {
    arclabels <- list()
    for (l in 1:len) {
      for (i in 1:(nnodes*nnodes)) {
        if (l==1) arclabels[[i]] <- sprintf("%.3f",lis[[l]][i])
        else arclabels[[i]][l] <- sprintf("%.3f",lis[[l]][i])
      }
    }
  }

  # construct a format string for nodes
  nodeformat <- ' %d [pos = "%f,%f", label = "%s", shape = %s, style = "%s",
    pin = true];'

  # construct a format string for edges
  if (usearclabels) {
    attributes <- 'style = %s, color = %s, label = "%s", dir = %s'
  }
  else {
    attributes <- 'style = %s, color = %s, dir = %s'
  }
  edgeformat = paste(sep='','  %d -> %d [', attributes, '];')

  # write the beginning of the file
  write('digraph G {', file=filename, append=FALSE, sep="")
  write('  center = 1;', file=filename, append=TRUE, sep="")
  write('  size = \"10, 10\";', file=filename, append=TRUE, sep="")
  write('  splines = true;', file=filename, append=TRUE, sep="")
  write('  sep = 1;', file=filename, append=TRUE, sep="")

  # write the nodes (observed and hidden)
  for (l in 1:(nrepeat+nlags)) {
    for (node in 1:nnodes) {
      thisnode <- sprintf(nodeformat, node+nnodes*(l-1), 100+200*(l-1),
                  100*(nnodes-node), nodelabels[[length(nodelabels)+1-
                  node-nnodes*(l-1)]],
                  nodeshapes[[node]], 'solid')
      write(thisnode, file=filename, append=TRUE, sep="")
    }
  }

  # write the nodes (closure)
  height <- seq(0,100*(nnodes-1),length.out=nclo)
  nl_closure <- NULL # nodelabes for closure variables
  for (i in 1:nclo) nl_closure <- c(nl_closure, paste(sep='','clo',i))
  for (node in 1:nclo) {
    thisnode <- sprintf(nodeformat, node+nnodes*(nrepeat+nlags), -100,
                height[nclo-node+1], nl_closure[node],
                nodeshapes[[node+nnodes]], 'solid')
    write(thisnode, file=filename, append=TRUE, sep="")
  }

  # write the edges (between observed and hidden varibables)
  threshold <- 0 # plot all edges 
             # 0.1 # only edges with connection absolute stronger than 0.1

  for (l in 1:len) {
    for (r in 1:(nrepeat+nlags-l+1)) { 
      for (node1 in 1:nnodes) {
        for (node2 in 1:nnodes) {
          if (abs(lis[[l]][node1,node2]) > threshold) {
            style <- 'solid'
            if(lis[[l]][node1,node2] > 0) color <- 'black' #'black' #'green'
            else color <- 'red'
            if (usearclabels) {
              thisarc <- sprintf(edgeformat, node2+(r-1)*nnodes,
                     node1+nnodes*(l-1)+(r-1)*nnodes, style, color,
                     arclabels[[(node2-1)*nnodes+node1]][l], 'forward')
            }
            else {
              thisarc <- sprintf(edgeformat, node2+(r-1)*nnodes,
                     node1+nnodes*(l-1)+(r-1)*nnodes, style, color, 'forward')
             }
             write(thisarc, file=filename, append=TRUE, sep="")
          } # end if (abs(...)>threshold)
        } # end node2
      } # end node1
    } # end r
  } # end l

  # write the edges (from the closure variables to other variables)
  for (i in 1:nclo) { # for each closure node have to write two edges
    color <- 'black' #'black' #'green'

    if (usearclabels) {
      thisarc <- sprintf(edgeformat, i+nnodes*(nrepeat+nlags),
              min(which(closure[,i]==1)), style, color,
              arclabels[[(node2-1)*nnodes+node1]][l], 'forward')
    }
    else {
      thisarc <- sprintf(edgeformat, i+nnodes*(nrepeat+nlags),
             min(which(closure[,i]==1)) , style, color, 'forward')
    }
    write(thisarc, file=filename, append=TRUE, sep="")

    if (usearclabels) {
      thisarc <- sprintf(edgeformat, i+nnodes*(nrepeat+nlags),
              max(which(closure[,i]==1)), style, color,
              arclabels[[(node2-1)*nnodes+node1]][l], 'forward')
    }
    else {
      thisarc <- sprintf(edgeformat, i+nnodes*(nrepeat+nlags),
             max(which(closure[,i]==1)) , style, color, 'forward')
    }
    write(thisarc, file=filename, append=TRUE, sep="")

  } # end i

  write('}', file=filename, append=TRUE, sep="")

}

