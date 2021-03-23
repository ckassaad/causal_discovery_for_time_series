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
    setwd("R_packages/VARLiNGAM/")
    # setwd('/home/kassaad/Documents/Codes/Survey_causal_discovery_time_series/scripts_R/R_packages/VARLiNGAM')
    source("VARLiNGAM.R")
    source("VAR_estim.R")
    source("ols_est.R")
    source("Gauss_Stats.R")
    source("Gauss_Tests.R")
    source("tsdata2canonicalform.R")
    # setwd('./lingam/code')
    source("./lingam/code/lingam.R")
    source("./lingam/code/estimate.R")
    source("./lingam/code/nzdiagbruteforce.R")
    source("./lingam/code/all.perm.R")
    source("./lingam/code/nzdiagscore.R")
    source("./lingam/code/iperm.R")
    source("./lingam/code/sltbruteforce.R")
    source("./lingam/code/sltscore.R")
    source("./lingam/code/prune.R")
    source("./lingam/code/tridecomp.R")
    source("./lingam/code/sqrtm.R")
}
start_up()

data_processed <- tsdata2canonicalform(data, nlags=n_lags)
result <- VARLiNGAM(data_processed, ntests=FALSE)
# result <- res$Mhat
result <- result$Bhat[[1]]

result_new <- data.frame(matrix(0, ncol = ncol(result), nrow = nrow(result)))
colnames(result_new) <- colnames(result)
rownames(result_new) <- colnames(result)
for (j1 in 1:nrow(result)){
    for (j2 in 1:nrow(result)){
      if (result[j1,j2] > alpha){
        result_new[j1,j2] <- 2
        result_new[j2,j1] <- 1
      }
    }
}

setwd(path)
#write.csv(result,"./results/result.csv", row.names = FALSE)
write.table(result_new, "./results/result.csv", col.names = colnames(data), row.names = colnames(data), sep = ",")
