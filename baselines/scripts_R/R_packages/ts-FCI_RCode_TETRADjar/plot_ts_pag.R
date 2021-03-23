# This function plots a time series pag
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



plot_ts_pag <- function(P, nobs=NULL, nodelabels=NULL,
    filename='tempfile') {

  # P ... a PAG consisting of all variables in the time serise (i.e. also the
  #       lagged variables). This is one big matrix, as returned by FCI.
  # nobs ... number of observed variables in 1 time slice
  
  if (is.null(nobs)) return("Please enter number of observed variables.")

  # generate the dot file
  dotfilename <- paste(sep='',filename,'.dot')
  psfilename <- paste(sep='',filename,'.ps')
 

  # write nodes and arrows to dot-file
  tspag2dot(P, nobs, dotfilename, nodelabels)


  # make a ps file using graphviz
  useprog <- 'neato' #-s -n2
  progcall <- paste(sep='',useprog,' -Tps -s72 ',dotfilename,' -o ',psfilename)
  system(progcall)

  # show the ps file using gv
  gvcall <- paste(sep='','gv ',psfilename)
  system(gvcall,wait=FALSE)

}


tspag2dot <- function(P, nobs, filename, nodelabels) {

  ntotal <- dim(P)[1]

  # if no nodelabels are given, create them as numbers
  if (is.null(nodelabels)){
    nodelabels <- list()
    for (i in 1:ntotal) {
      nodelabels[[i]] <- sprintf("%d",i)
    }
  }

  # construct a format string for nodes
  nodeformat <- ' %d [pos = "%f,%f", label = "%s", style = "%s", pin = true];'

  # construct a format string for edges
  attributes <- 'style = %s, color = %s, arrowtail = %s, arrowhead = %s'
  edgeformat = paste(sep='','  %d -> %d [', attributes, '];')

  # if filename is empty then write to a tempfile
  if (filename == "") {
    filename <- "tempfile.dot"
  }
  
  # write the beginning of the file
  write('digraph G {', file=filename, append=FALSE, sep="")
  write('  center = 1;', file=filename, append=TRUE, sep="")
  write('  size = \"10, 10\";', file=filename, append=TRUE, sep="")
  write('  splines = true;', file=filename, append=TRUE, sep="")
  write('  sep = 1;', file=filename, append=TRUE, sep="")

  # write the nodes
  cnt <- 0
  for (node in 1:ntotal) {
    if (node > 1 & (node-1)%%nobs == 0) cnt <- cnt + 1
    thisnode <- sprintf(nodeformat, node, 100+200*cnt,
                150*(ntotal-(node-1)%%nobs), nodelabels[[node]], 'solid')
    write(thisnode, file=filename, append=TRUE, sep="")
  }

  # write the edges
  arrows <- c("none", "normal", "odot") # styles of arrowheads and -tails
  Plow <- array(0, dim=dim(P))
  Plow[lower.tri(P)] <- P[lower.tri(P)] # each edge has two entries, only look
                                        # at each edge once
  for (node1 in 1:ntotal) {
    arcs <- which(Plow[node1,]!=0)
      for (node2 in arcs) {
        style <- 'solid'
        color <- 'black'

      thisarc <- sprintf(edgeformat, node2, node1, style, color,
                 arrows[P[node1,node2]], arrows[P[node2,node1]] )

        write(thisarc, file=filename, append=TRUE, sep="")
    } # end for node2
  } # end for node1

  write('}', file=filename, append=TRUE, sep="")

}
