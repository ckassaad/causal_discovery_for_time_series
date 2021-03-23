#'@method plot gimmep
#'@export
plot.gimmep <- function(x, file = NULL, subgroup = NULL, ...){
  if (is.null(file) & is.null(subgroup)) {
    cat("Please specify a file id for individual plots. Otherwise, summary plot is presented.")
    plot(x$group_plot_paths)
    invisible(x$group_plot_paths)
  } else if (is.null(subgroup)){
    a <- x$plots[[file]]
    plot(a)
    invisible(a)
  } else if (is.null(file)){
    ## insert code here to grab subgroup plot
    a <- x$sub_plots_paths[[subgroup]]
    if (!is.list(a)){
      cat("Subgroup", subgroup, "contains one individual. No subgroup plot provided.") 
    } else {
      plot(a)
      invisible(a)
    }
  }
}

#'@method print gimmep
#'@export
print.gimmep <- function(x, file = NULL, subgroup = NULL, 
                         mean = FALSE, estimates = FALSE, fitMeasures = FALSE, ...){
  if (!is.null(file)){
    if (estimates == TRUE){
      ind <- x$path_se_est[x$path_se_est$file %in% file, ]
      ind <- ind[, !names(ind) %in% c("subgroup")]
      ind[ ,4:6] <- round(ind[ ,4:6], digits = 3)
      cat("Coefficients for", file, "\n")
      print(ind, row.names = F)
      invisible(ind)
    } else if (estimates == FALSE){
      ind <- x$path_est_mats[[file]]
      colnames(ind) <- x$varnames
      # changed to reflect possibility of different # of lagged vars 
      # rownames(ind) <- x$varnames[(x$n_rois +1):(x$n_rois*2)]
      rownames(ind) <- x$varnames[(x$n_lagged+1):(x$n_vars_total)]
      ind <- round(ind, digits = 2)
      #ind_lag <- ind[ , 1:x$n_lagged]
      ind_lag <- ind[ , 1:x$n_lagged]
      #ind_con <- ind[ , (x$n_lagged+1):(x$n_rois*2)]
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
      allfit <- x$fit[, !names(x$fit) %in% c("modularity")]
      ind    <- allfit[allfit$file %in% file, ]
      cat("Fit for file", file, "\n")
      print.data.frame(ind, row.names = F)
    }
  }
  if (!is.null(subgroup)){
    if (mean == TRUE){
      # subidx   <- x$fit[,c("file", "subgroup")]
      # STL I'm not sure if this is what you meant: 
      subidx   <- x$fit[,c("file", "sub_membership")]
      files    <- subidx[subidx$sub_membership == subgroup, ]$file
      subfiles <- x$path_est_mats[files]
      if (length(files) == 1){
      cat("Subgroup", subgroup, "contains one individual. No average subgroup matrix provided.") 
      } else {
      s        <- apply(simplify2array(subfiles), 1:2, mean, na.rm = TRUE)
      colnames(s) <- x$varnames
      #rownames(s) <- x$varnames[(x$n_rois +1):(x$n_rois*2)]
      rownames(s) <- x$varnames[(x$n_lagged +1):(x$n_vars_total)]
      s     <- round(s, digits = 2)
      #s     <- s[(x$n_rois+1):(x$n_rois*2), ] #kmg did not comment this out 5/3/2018
      # s_lag <- s[ , 1:x$n_rois] #commented this out kmg
      s_lag <- s[ , 1:x$n_lagged]
      s_con <- s[ , (x$n_lagged+1):(x$n_vars_total)]
      cat("\n")
      cat("Lagged Average Matrix for Subgroup", subgroup, "\n")
      print(s_lag)
      cat("\n")
      cat("Contemporaneous Average Matrix for Subgroup", subgroup, "\n")
      print(s_con)
      invisible(s)
      }
    } else if (mean == FALSE & estimates == TRUE){
      # INSERT ESTIMATES FOR THAT SUBGROUP
      subest <- x$path_se_est
      subest <- subest[subest$subgroup == subgroup, ]
      subest[ ,4:6] <- round(subest[ ,4:6], digits = 3)
      cat("Coefficients for individuals in subgroup", subgroup, "\n")
      print(subest, row.names = F)
      invisible(subest)
    } else if (mean == FALSE & estimates == FALSE){
        sub <- x$path_counts_sub[[subgroup]]
        if (is.null(sub)){
          cat("Subgroup", subgroup, "contains one individual. No subgroup matrix provided.") 
        } else {
          #sub_lag <- sub[ , 1:x$n_rois]
          #sub_con <- sub[ , (x$n_rois+1):(x$n_rois*2)]
          sub_lag <- sub[ , 1:x$n_lagged]
          sub_con <- sub[ , (x$n_lagged+1):(x$n_vars_total)]
          cat("\n")
          cat("Lagged Count Matrix for Subgroup", subgroup, "\n")
          print(sub_lag)
          cat("\n")
          cat("Contemporaneous Count Matrix for Subgroup", subgroup, "\n")
          print(sub_con)
          invisible(sub)
        }
    }
    if (fitMeasures == TRUE){
      # RETURN FITMEASURES FOR MEMBERS OF THAT SUBGROUP
      allfit <- x$fit[, !names(x$fit) %in% c("modularity")]
      allfit[allfit$subgroup == subgroup, ]
      cat("Fit for individuals in subgroup", subgroup, "\n")
      print.data.frame(allfit, row.names = F)
    }
  }
  if (is.null(file) & is.null(subgroup)) {
    if (mean == FALSE & estimates == FALSE){
      cat("Please specify a file id for individual coefficient matrix. ", "\n", 
          "Otherwise, a summary count matrix is presented below.", "\n")
      all <- x$path_counts
      #all_lag <- all[ , 1:x$n_rois]
      all_lag <- all[ , 1:x$n_lagged]
      #all_con <- all[ , (x$n_rois+1):(x$n_rois*2)]
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
      #rownames(all2) <- x$varnames[(x$n_rois +1):(x$n_rois*2)]
      rownames(all2) <- x$varnames[(x$n_lagged +1):(x$n_vars_total)]
      all2 <- round(all2, digits = 2)
      #all2_lag <- all2[ , 1:x$n_rois]
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
      allfit <- x$fit[, !names(x$fit) %in% c("modularity")]
      cat("Fit for sample", "\n")
      print.data.frame(allfit, row.names = F)
      invisible(allfit)
    }
  }
}
  # else if (length(file) > 1){
  #   ind <- x$a[file]
  #   for (i in 1:length(file)){
  #     w <- ind[[i]]
  #     colnames(w) <- x$b
  #     rownames(w) <- x$b
  #     w <- round(w, digits = 2)
  #     w_lag <- w[(x$c+1):(x$c*2), 1:x$c]
  #     w_con <- w[(x$c+1):(x$c*2), (x$c+1):(x$c*2)]
  #     cat("\n")
  #     cat("Lagged Coefficient Matrix for", file[i], "\n")
  #     print(w_lag)
  #     cat("\n")
  #     cat("Contemporaneous Coefficient Matrix for", file[i], "\n")
  #     print(w_con)
  #   }
  #   invisible(w)
  # } 


print.gimme <- function(x, y, z){
  writeLines("gimme finished running normally")
  if (!is.null(z$out)) writeLines(paste("output is stored in", z$out))
  if (y == TRUE) {
    writeLines(paste("Number of subgroups =", x$n_subgroups))
    writeLines(paste("Modularity =", round(x$modularity, digits = 5)))
  }
}

