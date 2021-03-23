#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
name_dataframe = args[1]
data <- read.csv(file = name_dataframe)
name_alpha = args[2]
alpha <- scan(name_alpha, what=double())
name_n_lags = args[3]
n_lags <- scan(name_n_lags, what=double())
path <- args[4]
print(path)
print("#######################")
setwd(path)

start_up <- function() {
  setwd("R_packages/ts-FCI_RCode_TETRADjar/")
  source('dconnected.R')
  source('genData.R')
  source('main_tetrad_fci.R')
  source('plot_timeseries.R')
  source('Plotting_Commands_Barplots.R')
  source('plot_ts_pag.R')
  source('realData_tsfci.R')
  source('scores.R')
  source('Simulation_Commands.R')
  source('Simulations_data_cont.R')
  source('Simulations_data_disc.R')
  source('Simulations_graph.R')
  source('Tetrad_R_interact.R')
  source('ts_functions.R')
}
start_up()


result <- realData_tsfci(data=data, sig=alpha, nrep=n_lags, inclIE=FALSE, alg="tscfci", datatype="continuous", makeplot=FALSE)
temporal_names = c()
for (i in 1:n_lags){
  for (name in colnames(data)){
    temporal_names <- c(temporal_names, paste(name,i-1, sep = "_"))
  }
}
colnames(result) = temporal_names
setwd(path)
# %write.csv(result, file="./results/result.csv", row.names = FALSE)
write.table(result, "./results/result.csv", col.names = temporal_names, row.names = temporal_names, sep = ",")




