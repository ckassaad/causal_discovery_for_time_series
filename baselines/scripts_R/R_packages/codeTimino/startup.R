# Copyright (c) 2013-2013 Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms.
setwd('/home/kassaad/Documents/Codes/R - codes/onlineCodeTimino/codeTimino')

library(gam)
library(kernlab)
library(gptk)
# library(gptk)

source("granger_causality.R")
source("timino_causality.R")
source("./util/hammingDistance.R")
source("./util/indtestAll.R")
source("./util/indtestHsic.R")
source("./util/indtestPcor.R")
source("./util/TSindtest.R")
source("./util/fitting_ts.R")

