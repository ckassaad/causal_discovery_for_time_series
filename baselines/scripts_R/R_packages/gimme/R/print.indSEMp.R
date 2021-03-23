#'@method plot indSEMp
#'@export
plot.indSEMp <- function(x, file = NULL, ...){
  if (is.null(file)) {
    cat("Please specify a file id for individual plots. Otherwise, summary plot is presented.")
    plot(x$group_plot)
    invisible(x$group_plot)
  } else if (!is.null(file)){
    a <- x$plots[[file]]
    plot(a)
    invisible(a)
  }
}


#'@method print indSEMp
#'@export
print.indSEMp <- function(x, file = NULL, mean = FALSE, estimates = FALSE, fitMeasures = FALSE, ...){
  if (!is.null(file)){
    if (estimates == TRUE){
      ind <- x$path_est_mats[x$path_est_mats$file %in% file, ]
      ind[ ,4:6] <- round(ind[ ,4:6], digits = 3)
      cat("Coefficients for", file, "\n")
      print(ind, row.names = F)
      invisible(ind)
    } else if (estimates == FALSE){
      ind <- x$path_est_mats[[file]]
      colnames(ind) <- x$varnames
      rownames(ind) <- x$varnames[(x$n_lagged+1):(x$n_vars_total)]
      ind <- round(ind, digits = 2)
      ind_lag <- ind[ , 1:x$n_lagged]
      ind_con <- ind[ , (x$n_lagged+1):(x$n_vars_total)]
      cat("\n")
      cat("Lagged Matrix for", file, "\n")
      print(ind_lag)
      cat("\n")
      cat("Contemporaneous Matrix for", file, "\n")
      print(ind_con)
      invisible(ind)
    }
    if (fitMeasures == TRUE){
      allfit <- x$fit
      ind    <- allfit[allfit$file %in% file, ]
      cat("Fit for file", file, "\n")
      print.data.frame(ind, row.names = F)
    }
  }
  if (is.null(file)) {
    if (mean == FALSE & estimates == FALSE){
      cat("Please specify a file id for individual coefficient matrix. ", "\n", 
          "Otherwise, a summary count matrix is presented below.", "\n")
      all <- x$path_counts
      all_lag <- all[ , 1:x$n_lagged]
      all_con <- all[ , (x$n_lagged+1):(x$n_vars_total)]
      cat("\n")
      cat("Lagged Count Matrix for Sample", "\n")
      print(all_lag)
      cat("\n")
      cat("Contemporaneous Count Matrix for Sample", "\n")
      print(all_con)
      invisible(all)
    } else if (mean == FALSE & estimates == TRUE){
      print(x$path_se_est) 
      invisible(x$path_se_est)
    }
    if (mean == TRUE){
      cat("Please specify a file id for individual coefficient matrix. ", "\n", 
          "Otherwise, a summary average matrix is presented below.", "\n")
      all2 <- apply(simplify2array(x$path_est_mats), 1:2, mean, na.rm = TRUE)
      colnames(all2) <- x$varnames
      rownames(all2) <- x$varnames[(x$n_lagged+1):(x$n_vars_total)]
      all2 <- round(all2, digits = 2)
      all2_lag <- all2[ , 1:x$n_lagged]
      all2_con <- all2[ , (x$n_lagged+1):(x$n_vars_total)]
      cat("\n")
      cat("Lagged Average Matrix for Sample", "\n")
      print(all2_lag)
      cat("\n")
      cat("Contemporaneous Average Matrix for Sample", "\n")
      print(all2_con)
    }
    if (fitMeasures == TRUE){
      allfit <- x$fit
      cat("Fit for sample", "\n")
      print.data.frame(allfit, row.names = F)
      invisible(allfit)
    }
  }
}



