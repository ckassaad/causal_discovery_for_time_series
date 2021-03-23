#' Grabs final coefficients for each individual.
#' @param dat A list containing information created in setup().
#' @param grp A list containing group-level information. NULL in aggSEM and
#' indSEM.
#' @param ind A list containing individual- and (potentially) subgroup-level
#' information.
#' @param k The counter indicating the individual.
#' @return Individual-level information on fit, coefficients, and plots.
#' @keywords internal
get.params <- function(dat, grp, ind, k){
  
  op  = NULL # appease CRAN check
  ind_plot = NA
  
  if (!dat$agg){
    fit <- fit.model(syntax    = c(dat$syntax, 
                                   grp$group_paths, 
                                   ind$sub_paths[[k]], 
                                   ind$ind_paths[[k]]), 
                     data_file = dat$ts_list[[k]])
  } else {
    data_file <- do.call("rbind", dat$ts_list)
    fit        <- fit.model(syntax    = c(dat$syntax, 
                                          grp$group_paths,
                                          ind$sub_paths[[k]], 
                                          ind$ind_paths[[k]]), 
                            data_file = data_file)
  }
  
  error   <- any(grepl("error", class(fit)))
  
  if (!error) {
    converge <- lavInspect(fit, "converged")
    zero_se  <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
  } else {
    converge <- FALSE
    zero_se  <- TRUE
  }
  
  # if no convergence, roll back one path at individual level, try again 
  if (!converge | zero_se){
    status <- "nonconvergence"
    if (length(ind$ind_paths[[k]]!= 0)){
      ind$ind_paths[[k]] <- ind$ind_paths[[k]][-length(ind$ind_paths[[k]])]
      if (!dat$agg){
        fit <- fit.model(syntax    = c(dat$syntax, 
                                       grp$group_paths, 
                                       ind$sub_paths[[k]], 
                                       ind$ind_paths[[k]]), 
                         data_file = dat$ts_list[[k]])
      } else {
        data_file <- do.call("rbind", dat$ts_list)
        fit        <- fit.model(syntax    = c(dat$syntax, 
                                              grp$group_paths,
                                              ind$sub_paths[[k]], 
                                              ind$ind_paths[[k]]), 
                                data_file = data_file)
      }
    }
    
    error   <- any(grepl("error", class(fit)))
    
    if (!error){
      converge  <- lavInspect(fit, "converged")
      
      ind_coefs0 <- standardizedSolution(fit)
      ind_coefs_idx <- paste0(ind_coefs0$lhs,ind_coefs0$op,ind_coefs0$rhs)
      ind_coefs <- ind_coefs0[ind_coefs0$op == "~" |
                                ind_coefs_idx %in% c(dat$candidate_paths, dat$candidate_corr),]
      #ind_coefs <- ind_coefs0[ind_coefs_idx %in% elig_paths,]
      #commented out by lan 4.11.2019
      #ind_coefs <- subset(standardizedSolution(fit), op == "~") # if betas = 0, no SEs
      if (length(ind_coefs[,1]) > 0){
        zero_se   <- sum(lavInspect(fit, "se")$beta, na.rm = TRUE) == 0
      } else {
        zero_se <- FALSE
      }
      if (converge){
        status <- "last known convergence"
      }
    } else {
      converge <- FALSE
      zero_se  <- TRUE
    }
  }
  
  if (converge & !zero_se){#& (ind$n_ind_paths[k] >0) ){
    status   <- "converged normally"
    
    ind_fit    <- fitMeasures(fit, c("chisq", "df", "npar", "pvalue", "rmsea", 
                                     "srmr", "nnfi", "cfi", "bic", "aic", "logl"))
    ind_fit    <- round(ind_fit, digits = 4)
    ind_fit[2] <- round(ind_fit[2], digits = 0)
    
    ind_vcov_full <- lavInspect(fit, "vcov.std.all")
    keep          <- rownames(ind_vcov_full) %in% dat$candidate_paths
    ind_vcov      <- ind_vcov_full[keep, keep]
    
    
    ind_coefs0 <- standardizedSolution(fit)
    ind_coefs_idx <- paste0(ind_coefs0$lhs,ind_coefs0$op,ind_coefs0$rhs)
    ind_coefs <- ind_coefs0[ind_coefs0$op == "~" |
                              ind_coefs_idx %in% c(dat$candidate_paths, dat$candidate_corr),]
    #ind_coefs <- subset(standardizedSolution(fit), op == "~")
    
    # if (length(ind_coefs[,1]) > 0){ # stl comment out 11.20.17
    ind_betas <- round(lavInspect(fit, "std")$beta, digits = 4)
    ind_ses   <- round(lavInspect(fit, "se")$beta, digits = 4)

    #added to ensure correct ordering in matrices
    ind_betas <- ind_betas[dat$varLabels$orig,]
    ind_betas <- ind_betas[,dat$varLabels$coln]
    
    ind_ses <- ind_ses[dat$varLabels$orig,]
    ind_ses <- ind_ses[,dat$varLabels$coln]
    
    # zf added 2019-01-23
    ind_psi <- round(lavInspect(fit, "std")$psi, digits = 4)
    ind_psi_unstd <- round(lavInspect(fit, "estimates")$psi, digits = 4)
    
    #rownames(ind_betas) <- rownames(ind_ses) <- dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)]
    #colnames(ind_betas) <- colnames(ind_ses) <- dat$varnames
    #   } # stl comment out 11.20.17 
    
    if (dat$agg & !is.null(dat$out)){
      
      write.csv(ind_betas, file.path(dat$out, "allBetas.csv"), 
                row.names = TRUE)
      
      # write.csv(ind_vcov_full, file.path(dat$out, "allvcov.csv"), 
      #           row.names = TRUE)
      
      write.csv(ind_ses, file.path(dat$out, "allStdErrors.csv"), 
                row.names = TRUE)
      
      # zf added 2019-01-23
      write.csv(ind_psi, file.path(dat$out, "allPsi.csv"),row.names = TRUE)
      write.csv(ind_psi_unstd, file.path(dat$out, "allPsiUnstd.csv"),row.names = TRUE)
      
    } else if (!dat$agg & !is.null(dat$out)) { # & ind$n_ind_paths[k]>0)
      write.csv(ind_betas, file.path(dat$ind_dir, 
                                     paste0(dat$file_order[k,2], 
                                            "Betas.csv")), row.names = TRUE)
      
      # write.csv(ind_vcov_full, file.path(dat$ind_dir, 
      #                                paste0(dat$file_order[k,2], 
      #                                       "vcov.csv")), row.names = TRUE)
      # zf added 2019-01-23
      write.csv(ind_psi, file.path(dat$ind_dir, 
                                     paste0(dat$file_order[k,2], 
                                            "Psi.csv")), row.names = TRUE)
      write.csv(ind_psi_unstd, file.path(dat$ind_dir, 
                                     paste0(dat$file_order[k,2], 
                                            "PsiUnstd.csv")), row.names = TRUE)
      write.csv(ind_ses, file.path(dat$ind_dir,
                                   paste0(dat$file_order[k,2], 
                                          "StdErrors.csv")), row.names = TRUE)
    }
    
    ind_plot  <- NA
    if (dat$plot){
      ind_betas_t <- t(ind_betas)
      lagged      <- ind_betas_t[1:dat$n_lagged, ]
      contemp     <- ind_betas_t[(dat$n_lagged+1):(dat$n_vars_total), ]
      plot_vals   <- rbind(w2e(lagged), w2e(contemp))
      is_lagged   <- c(rep(TRUE, sum(lagged != 0)), 
                       rep(FALSE, sum(contemp != 0)))
      
      plot_file   <- ifelse(dat$agg, 
                            file.path(dat$out, "summaryPathsPlot.pdf"),
                            file.path(dat$ind_dir, 
                                      paste0(dat$file_order[k,2], "Plot.pdf")))
      
      ind_plot <- tryCatch(qgraph(plot_vals,
                                  layout       = "circle",
                                  lty          = ifelse(is_lagged, 2, 1),
                                  edge.labels  = FALSE,
                                  curve        = FALSE,
                                  parallelEdge = TRUE,
                                  fade         = FALSE,
                                  posCol       = "red",
                                  negCol       = "blue",
                                  labels       = 
                                    dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                                  label.cex    = 2,
                                  DoNotPlot    = TRUE), 
                           error = function(e) e)
      
      if (!is.null(dat$out) & !"error" %in% class(ind_plot)){
        pdf(plot_file)
        plot(ind_plot)
        dev.off()
      }
    }
  } 
  
  # commented out on 11.20.17 by stl 
  # if (ind$n_ind_paths[k] ==0 & converge) {
  #   status     <- "no paths added"
  #   ind_fit    <- fitMeasures(fit, c("chisq", "df", "npar", "pvalue", "rmsea", 
  #                                    "srmr", "nnfi", "cfi", "bic", "aic", "logl"))
  #   ind_fit    <- round(ind_fit, digits = 4)
  #   ind_fit[2] <- round(ind_fit[2], digits = 0)
  #   
  #   ind_vcov  <- lavInspect(fit, "vcov.std.all")
  #   keep      <- rownames(ind_vcov) %in% dat$candidate_paths
  #   ind_vcov  <- ind_vcov[keep, keep]
  #   
  #   ind_betas <- NULL
  #   ind_coefs <- subset(standardizedSolution(fit), op == "~")
  #   
  # } 
  
  if (!converge | zero_se){
    if (!converge) status <- "nonconvergence"
    if (zero_se)   status <- "computationally singular"
    ind_fit   <- rep(NA, 11)
    ind_coefs <- matrix(NA, nrow = 1, ncol = 9)
    colnames(ind_coefs) <- c("lhs", "op", "rhs", "est.std", "se", "z", "pvalue", "ci.lower", "ci.upper")
    ind_betas <- NA
    ind_vcov  <- NA
    ind_plot  <- NA
    ind_psi   <- NA
    ind_psi_unstd   <- NA
    ind_vcov_full <- NA
  }
  
  if (!dat$plot)
    ind_plot  <- NA
  
  res <- list("status"    = status, 
              "ind_fit"   = ind_fit, 
              "ind_coefs" = ind_coefs, 
              "ind_betas" = ind_betas, 
              "ind_psi"   = ind_psi, 
              "ind_psi_unstd" = ind_psi_unstd, 
              "ind_vcov"  = ind_vcov,
              "ind_vcov_full"  = ind_vcov_full,
              "ind_plot"  = ind_plot,
              "ind_syntax" = c(dat$syntax, grp$group_paths,ind$sub_paths[[k]], ind$ind_paths[[k]])
              )
  return(res)
}
