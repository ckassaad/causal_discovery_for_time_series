#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
name_dataframe = args[1]
data <- read.csv(file = name_dataframe)
data$time_index <- NULL
name_alpha = args[2]
alpha <- scan(name_alpha, what=double())
name_n_lags = args[3]
n_lags <- scan(name_n_lags, what=double())
path <- args[4]
print(path)
print("#######################")
setwd(path)

start_up <- function() {
    library(gam)
    library(kernlab)
    library(gptk)
    source("./R_packages/codeTimino/timino_causality.R")
    source("./R_packages/codeTimino/util/hammingDistance.R")
    source("./R_packages/codeTimino/util/indtestAll.R")
    source("./R_packages/codeTimino/util/indtestHsic.R")
    source("./R_packages/codeTimino/util/indtestPcor.R")
    source("./R_packages/codeTimino/util/TSindtest.R")
    source("./R_packages/codeTimino/util/fitting_ts.R")
}
start_up()

result <- timino_dag(data, alpha = alpha, max_lag = n_lags, model = traints_linear, indtest = indtestts_crosscov, output = TRUE)

result[is.na(result)] <- 3

for (j1 in 1:nrow(result)){
    for (j2 in 1:nrow(result)){
      if (result[j1,j2] == 1){
        result[j1,j2] <- 2

      }
    }
}

for (j1 in 1:nrow(result)){
    for (j2 in 1:nrow(result)){
      if (result[j1,j2] == 2){
        if (result[j2,j1] == 0){
            result[j2,j1] <- 1
        }
      }
      if (j1 == j2){
          result[j1,j2] <- 1
      }
    }
}

for (j1 in 1:nrow(result)){
    for (j2 in 1:nrow(result)){
      if (result[j1,j2] == 3){
        result[j1,j2] <- 2
        result[j2,j1] <- 2
      }
    }
}


#write.csv(result, file="./results/result.csv", row.names = FALSE, col.names = FALSE)
write.table(result, "./results/result.csv", col.names = colnames(data), row.names = colnames(data), sep = ",")

