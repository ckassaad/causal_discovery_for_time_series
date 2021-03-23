#'@method plot aggSEMp
#'@export
plot.aggSEMp <- function(x, ...){
    plot(x$plot[[1L]])
    invisible(x$plot[[1L]])
}

#'@method print aggSEMp
#'@export
print.aggSEMp <- function(x, estimates = FALSE, fitMeasures = FALSE, ...){
    if (estimates == TRUE){
      cat("Coefficients for final model", "\n")
      print(x$path_se_est, row.names = F)
      invisible(x$path_se_est)
    } else if (estimates == FALSE){
      ind <- x$path_est_mat[[1]]
      colnames(ind) <- x$varnames
      rownames(ind) <- x$varnames[(x$n_lagged+1):(x$n_vars_total)]
      ind <- round(ind, digits = 2)
      ind_lag <- ind[ , 1:x$n_lagged]
      ind_con <- ind[ , (x$n_lagged+1):(x$n_vars_total)]
      cat("\n")
      cat("Lagged Matrix for all", "\n")
      print(ind_lag)
      cat("\n")
      cat("Contemporaneous Matrix for all", "\n")
      print(ind_con)
      invisible(ind)
    }
    if (fitMeasures == TRUE){
      cat("Fit for all", "\n")
      print.data.frame(x$fit, row.names = F)
    }
  }
  