# To call ts(C)FCI on a set of time series data, calls TETRAD
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



realData_tsfci <- function(data=NULL, path=NULL, sig=0.01, nrep=4, inclIE=FALSE,
  alg="tscfci", datatype="continuous", makeplot=FALSE) {

  # either give data in a matrix or a path to a data file (in txt-format) from
  # where the data should be loaded:
  # data ... time series data, each column is one variable and each row is one
  #          observation (data should thus have much more rows than columns)
  # path ... string of the path to a txt file containing the data (same format
  #          as above)

  # sig ... significance level for independence test in ts(c)fci
  # nrep ... lenght of the unrolled time series graph (how often the _nodes_
  #          are included), nrep = tau + 1 (window length in paper)
  # inclIE ... boolean, if FALSE no instantaneous effects are included (that
  #          means that all intantaneous edges are oriented as double headed
  #          arrows, ie. are due to confounders;, if TRUE, instantaneous 
  #          effects are inlcuded
  # alg ... "tsfci" of "tscfci" (it's recommended to use "tscfci")
  # datatype ... "continuous" or "discrete" depending on the data
  # makeplot ... boolean, if TRUE the output PAG of ts(C)FCI is plotted, needs
  #          the graphviz program

  if (is.null(data) & is.null(path)) {
    return("Input either data-matrix or path to txt data-file")
  }

  # load Data
  if (is.null(data)) {
    print("Loading data.")
    data <- read.table(path,header=TRUE)
  }

  nobs <- ncol(data)

  # write data to file
  print("Write data to file.")
  writeDataToFile(data, nrep)

  # write background knowledge to file
  print("Write knowledge to file.")
  writeKnowledgeToFile(nobs, nrep)

  # run ts(c)fci from Tetrad and save output PAG to "resultFCI.txt"
  print("Call Tetrad.")
  callTetradUniversal("tsdata", algName=alg, dataType=datatype, sig=sig,
    inclInstEffect=inclIE)

  # read in results from Fci
  print("Read in ts(C)FCI output.")
  pag <- readFCIoutput()

  if (makeplot) {
    # plot the output PAG
    print("Plot PAG.")
    plot_ts_pag(pag,nobs)
  }

  pag

}
