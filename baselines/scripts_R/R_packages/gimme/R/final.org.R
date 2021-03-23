#' Wrapup, create output files.
#' @param dat A list containing information created in setup().
#' @param grp A list containing group-level information. NULL in aggSEM and
#' indSEM.
#' @param ind A list containing individual- and (potentially) subgroup-level
#' information.
#' @param sub A list containing subgroup information.
#' @param sub_spec A list containing information specific to each subgroup.
#' @param store A list containing output from indiv.search().
#' @return Aggregated information, such as counts, levels, and plots.
#' @keywords internal
final.org <- function(dat, grp, sub, sub_spec, diagnos=FALSE, store){
  
  ind = store$ind
  sub_coefs  <- list()
  sub_summ   <- list()
  sub_plots  <- list()
  sub_plots_cov <- list()
  sub_counts <- list()
  sub_counts_cov <- list()
  
  param  = NULL # appease CRAN check
  est.std = NULL # appease CRAN check
  count  = NULL # appease CRAN check
  
  if (!dat$agg){
    
    coefs       <- do.call("rbind", store$coefs)
    
    if(length(coefs[,1])>0){
      coefs$id    <- rep(names(store$coefs), sapply(store$coefs, nrow))
      coefs$param <- paste0(coefs$lhs, coefs$op, coefs$rhs)
      coefs <- coefs[!coefs$param %in% dat$nonsense_paths,] # Removes non-sense paths that occur when ar = FALSE or mult_vars is not null from output 
      
      coefs$level[coefs$param %in% c(grp$group_paths, dat$syntax)] <- "group"
      coefs$level[coefs$param %in% unique(unlist(ind$ind_paths))]  <- "ind"
      coefs$color[coefs$level == "group"] <- "black"
      coefs$color[coefs$level == "ind"]   <- "gray50"
    }
    
    indiv_paths <- NULL
    samp_plot <- NULL
    samp_plot_cov <- NULL
    sample_counts <- NULL
    sample_counts_corr <- NULL
    # if (length(coefs[,1])>0){ # commented out stl 11.20.17
    if (dat$subgroup) {
      if (sub$n_subgroups != dat$n_subj){ # ensure everyone isn't in their own subgroup
        
        sub_paths_count <- table(unlist(
          lapply(sub_spec, FUN = function(x) c(x$sub_paths))))
        # if path exists in all subgroups, push it up to group level
        sub_to_group    <- names(
          sub_paths_count[sub_paths_count == sub$n_subgroups])
        
        
        for (s in 1:sub$n_subgroups){
          sub_s_mat_counts <- matrix(0, nrow = (dat$n_vars_total), 
                                     ncol = (dat$n_vars_total))
          sub_s_mat_counts_cov <- sub_s_mat_counts
          sub_s_mat_means  <- sub_s_mat_counts
          sub_s_mat_means_cov  <- sub_s_mat_counts
          sub_s_mat_colors <- matrix(NA, nrow = (dat$n_vars_total), 
                                     ncol = (dat$n_vars_total))
          sub_s_mat_colors_cov <- matrix(NA, nrow = (dat$n_vars_total), 
                                     ncol = (dat$n_vars_total))
          
          sub_s_coefs <- coefs[coefs$id %in% sub_spec[[s]]$sub_s_subjids, ]
          sub_s_coefs$level[sub_s_coefs$param %in% sub_spec[[s]]$sub_paths] <- "sub"
          sub_s_coefs$level[sub_s_coefs$param %in% sub_to_group] <- "group"
          sub_s_coefs$color[sub_s_coefs$level == "group"] <- "black"
          sub_s_coefs$color[sub_s_coefs$level == "sub"]   <- "green3"
          sub_s_coefs$color[sub_s_coefs$level == "ind"]   <- "gray50"
          
          ## march 2018 stl - fix to remove error caused where lhs and rhs 
          ## values are NA. there's no deeper trouble here - it was just due to an 
          ## rbind where individuals with no paths (e.g., entirely NA) were included
          ## in the full rbind, which led to variable names of "NA" 
          sub_s_coefs <- sub_s_coefs[!is.na(sub_s_coefs$lhs), ]
          sub_s_coefs <- sub_s_coefs[!is.na(sub_s_coefs$rhs), ]
          
          sub_s_summ <- transform(sub_s_coefs, 
                                  count = as.numeric(
                                    ave(param, param, FUN = length)),
                                  mean  = ave(est.std, param, FUN = sum)/sub_spec[[s]]$n_sub_subj)
          sub_s_summ <- subset(sub_s_summ, !duplicated(param))
          sub_s_summ$row <- match(sub_s_summ$lhs, dat$lvarnames)
          sub_s_summ$col <- match(sub_s_summ$rhs, dat$lvarnames)
          sub_s_summ$mem <- s
          
          regressions <- 
            sub_s_summ[which(sub_s_summ$op == "~"),]
          sub_s_mat_counts[cbind(regressions$row, regressions$col)] <- 
            as.numeric(as.character(regressions$count))
          sub_s_mat_counts <- sub_s_mat_counts[(dat$n_lagged+1):(dat$n_vars_total), ]
          colnames(sub_s_mat_counts) <- dat$varnames
          rownames(sub_s_mat_counts) <- dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)]
          
          sub_s_mat_means[cbind(regressions$row, regressions$col)]  <- regressions$mean
          sub_s_mat_colors[cbind(regressions$row, regressions$col)] <- regressions$color
          sub_s_mat_colors <- sub_s_mat_colors[(dat$n_lagged+1):(dat$n_vars_total), ]
          
          cov <- 
            sub_s_summ[(which(sub_s_summ$op == "~~")),]
          sub_s_mat_counts_cov[cbind(cov$row, cov$col)] <- 
            as.numeric(as.character(cov$count))
          sub_s_mat_counts_cov <- sub_s_mat_counts_cov[(dat$n_lagged+1):(dat$n_vars_total), ]
          colnames(sub_s_mat_counts_cov) <- dat$varnames
          rownames(sub_s_mat_counts_cov) <- dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)]
          
          sub_s_mat_means_cov[cbind(cov$row, cov$col)]  <- cov$mean
          sub_s_mat_colors_cov[cbind(cov$row, cov$col)] <- cov$color
          sub_s_mat_colors_cov <- sub_s_mat_colors_cov[(dat$n_lagged+1):(dat$n_vars_total), ]
          
          if (dat$plot & sub_spec[[s]]$n_sub_subj != 1){ #plot subgroup plot if >1 nodes in subgroup
            
            sub_s_counts <- t(sub_s_mat_counts/sub_spec[[s]]$n_sub_subj)
            lagged     <- sub_s_counts[1:(dat$n_lagged), ]

            contemp    <- sub_s_counts[(dat$n_lagged+1):(dat$n_vars_total), ]
            plot_vals  <- rbind(w2e(lagged), w2e(contemp))
            is_lagged  <- c(rep(TRUE, sum(lagged != 0)), rep(FALSE, sum(contemp != 0)))
            
            sub_colors <- t(sub_s_mat_colors)
            colors     <- c(sub_colors[1:(dat$n_lagged), ],
                            sub_colors[(dat$n_lagged+1):(dat$n_vars_total), ])
            colors     <- colors[!is.na(colors)]
            
            sub_plot <- tryCatch(qgraph(plot_vals,
                                        layout       = "circle",
                                        lty          = ifelse(is_lagged, 2, 1),
                                        edge.labels  = FALSE,
                                        edge.color   = colors,
                                        parallelEdge = TRUE,
                                        fade         = FALSE,
                                        labels       = 
                                          dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                                        label.cex    = 2,
                                        DoNotPlot    = TRUE), 
                                 error = function(e) e)
            
            sub_s_counts_cov <- t(sub_s_mat_counts_cov/sub_spec[[s]]$n_sub_subj)
            contemp_cov    <- sub_s_counts_cov[(dat$n_lagged+1):(dat$n_vars_total), ]
            plot_vals_cov  <- w2e(contemp_cov)
            sub_colors_cov <- t(sub_s_mat_colors_cov)
            #commented out by lan 2.10.2020
            #colors     <- c(sub_colors_cov[1:(dat$n_lagged), ],
            #sub_colors_cov[(dat$n_lagged+1):(dat$n_vars_total), ])
            colors    <- sub_colors_cov[(dat$n_lagged+1):(dat$n_vars_total), ]
            colors     <- colors[!is.na(colors)]
            sub_plot_cov <- tryCatch(qgraph(plot_vals_cov,
                                              layout       = "circle",
                                              edge.labels  = FALSE,
                                              edge.color   = colors,
                                              parallelEdge = TRUE,
                                              fade         = FALSE,
                                              arrows       = FALSE,
                                              labels       = 
                                                dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                                              label.cex    = 2,
                                              DoNotPlot    = TRUE), 
                                       error = function(e) e)
            
            if (!is.null(dat$out) & !"error" %in% class(sub_plot)){
              pdf(file.path(dat$subgroup_dir, 
                            paste0("subgroup", s, "Plot.pdf")))
              plot(sub_plot)
              dev.off()
              if(sum(sub_s_counts_cov)>0){
              pdf(file.path(dat$subgroup_dir, 
                            paste0("subgroup", s, "Plot_cov.pdf")))
              plot(sub_plot_cov)
              dev.off()}
            }
            
          } else {
            sub_plot         <- NULL
            sub_s_mat_counts <- NULL
          }
          
          if (sub_spec[[s]]$n_sub_subj != 1 & !is.null(dat$out)){
            write.csv(sub_s_mat_counts, 
                      file = file.path(dat$subgroup_dir, 
                                       paste0("subgroup", s, 
                                              "PathCountsMatrix.csv")), 
                      row.names = TRUE)
          }
          if (dat$plot & sub_spec[[s]]$n_sub_subj != 1){ ##add by lan 021220: store the sub_plot & sub_plot_cov to the plots when n>1
          sub_coefs[[s]] <- sub_s_coefs
          sub_summ[[s]]  <- sub_s_summ
          sub_plots[[s]] <- sub_plot
          sub_plots_cov[[s]] <- sub_plot_cov
          sub_counts[[s]] <- sub_s_mat_counts
          sub_counts_cov[[s]] <- sub_s_mat_counts_cov
          }
        }
        
        summ <- do.call("rbind", sub_summ)
        coefs <- do.call("rbind", sub_coefs)
        
      } else {
        sub_coefs <- NULL
        sub_plots <- NULL
        sub_paths <- NULL
        sub_plots_cov <- NULL
        sub_counts_cov <- NULL
        summ <- transform(coefs, count = as.numeric(
          ave(param, param, FUN = length)))
        summ <- subset(summ, !duplicated(param)) 
      }
    }
    else {
      sub_coefs <- NULL
      sub_plots <- NULL
      sub_paths <- NULL
      sub_plots_cov <- NULL
      sub_counts_cov <- NULL
      summ <- transform(coefs, count = as.numeric(
        ave(param, param, FUN = length)))
      summ <- subset(summ, !duplicated(param)) 
    }
    
    # combining and creating wide summaryPathCounts -------------------------- #
    summ$label <- ifelse(summ$level == "sub", 
                         paste0("subgroup", summ$mem),
                         summ$level)
    a <- aggregate(count ~ lhs + op + rhs + label, data = summ, sum)
    
    a <- a[order(-a$count, a$label),]
    a <- reshape(a, timevar = "label", idvar = c("lhs", "op", "rhs"), 
                 direction = "wide")
    a[is.na(a)] <- 0
    a$lhs <- recode.vars(a$lhs, dat$lvarnames, dat$varnames)
    a$rhs <- recode.vars(a$rhs, dat$lvarnames, dat$varnames)
    
    if (!is.null(dat$out)){
      write.csv(a, file.path(dat$out, "summaryPathCounts.csv"), 
                row.names = FALSE)
    }
    
    # end creating wide summaryPathCounts ------------------------------------ #
    
    b <- aggregate(count ~ lhs + op + rhs + color + label + param, data = summ, sum)
    b <- transform(b, xcount = ave(count, param, FUN = sum))
    # sorting by count and then dropping duplicated parameters
    # ensures that subgroup paths will be displayed as green
    # and individual paths that appear in another subgroup
    # will not cause subgroup paths to all display as individual
    # CA 10.5.18 created variable to order by label.  Some individual paths were 
    # being selected over subgroup paths in the duplicated function.
    
    b$labelnum[b$label=='group'] <- 1
    b$labelnum[b$label=='ind'] <- 3
    b$labelnum[is.na(b$labelnum)] <-2
    
    b <- b[order(b$labelnum), ]
    d <- b[!duplicated(b$param), c("lhs", "op", "rhs", "color", "xcount")] 
    
    c_direct <- d[which(d$op == "~"),]
    c_corr   <- d[which(d$op == "~~"),]
    
    c_direct$row <- match(c_direct$lhs, dat$lvarnames) - dat$n_lagged
    c_direct$col <- match(c_direct$rhs, dat$lvarnames)
    c_corr$row <- match(c_corr$lhs, dat$lvarnames) - dat$n_lagged
    c_corr$col <- match(c_corr$rhs, dat$lvarnames)
    
    sample_counts <- matrix(0, ncol = (dat$n_vars_total), nrow = (dat$n_vars_total - dat$n_lagged))
    sample_counts[cbind(c_direct$row, c_direct$col)] <- c_direct$xcount
    sample_counts_corr <- matrix(0, ncol = (dat$n_vars_total), nrow = (dat$n_vars_total - dat$n_lagged))
    sample_counts_corr[cbind(c_corr$row, c_corr$col)] <- c_corr$xcount
    colnames(sample_counts) <- dat$varnames
    rownames(sample_counts) <- dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)]
    colnames(sample_counts_corr) <- dat$varnames
    rownames(sample_counts_corr) <- dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)]
    
    if (dat$plot){
      
      sample_colors <- matrix(NA, ncol = (dat$n_vars_total), nrow = (dat$n_vars_total-dat$n_lagged))
      sample_colors[cbind(c_direct$row, c_direct$col)] <- c_direct$color
      sample_colors_corr <- matrix(NA, ncol = (dat$n_vars_total), nrow = (dat$n_vars_total-dat$n_lagged))
      sample_colors_corr[cbind(c_corr$row, c_corr$col)] <- c_corr$color
      
      sample_paths  <- t(sample_counts)/dat$n_subj
      sample_paths_corr <- t(sample_counts_corr)/dat$n_subj
      
      lagged     <- sample_paths[1:(dat$n_lagged), ]
     
      contemp    <- sample_paths[(dat$n_lagged+1):(dat$n_vars_total), ]
      plot_vals  <- rbind(w2e(lagged), w2e(contemp))
      is_lagged  <- c(rep(TRUE, sum(lagged != 0)),
                      rep(FALSE, sum(contemp != 0)))
      
      samp_colors <- t(sample_colors)
      colors      <- c(samp_colors[1:(dat$n_lagged), ],
                       samp_colors[(dat$n_lagged+1):(dat$n_vars_total), ])
      colors      <- colors[!is.na(colors)]
      
      samp_plot <- tryCatch(qgraph(plot_vals,
                                   layout       = "circle",
                                   lty          = ifelse(is_lagged, 2, 1),
                                   edge.labels  = FALSE,
                                   edge.color   = colors,
                                   parallelEdge = TRUE,
                                   fade         = FALSE,
                                   labels       = 
                                     dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                                   label.cex    = 2,
                                   DoNotPlot    = TRUE), 
                            error = function(e) e)
      samp_colors_corr <- t(sample_colors_corr)
      #commented out by lan 2.10.2020
      # colors_corr      <- c(samp_colors_corr[1:(dat$n_lagged), ],
      #                samp_colors_corr[(dat$n_lagged+1):(dat$n_vars_total), ])
      colors_corr     <- samp_colors_corr[(dat$n_lagged+1):(dat$n_vars_total)]
      colors_corr      <- colors_corr[!is.na(colors_corr)]
      
      if (sum(sample_paths_corr)>0){
        corr   <- sample_paths_corr[(dat$n_lagged+1):(dat$n_vars_total), ]
        plot_vals_corr  <- w2e(corr)
        samp_plot_cov <- tryCatch(qgraph(plot_vals_corr,
                                          layout       = "circle",
                                          edge.labels  = FALSE,
                                          edge.color   = colors_corr,
                                          parallelEdge = TRUE,
                                          fade         = FALSE,
                                          arrows       = FALSE,
                                          labels       = 
                                            dat$varnames[(dat$n_lagged+1):(dat$n_vars_total)],
                                          label.cex    = 2,
                                          DoNotPlot    = TRUE), 
                                   error = function(e) e)
      }
       
      
      if (!is.null(dat$out)){
        pdf(file.path(dat$out, "summaryPathsPlot.pdf"))
        plot(samp_plot)
        dev.off()
        if(sum(sample_paths_corr)>0){
        pdf(file.path(dat$out, "summaryCovPlot.pdf"))
        plot(samp_plot_cov)
        dev.off()}
      }
      
    } else {
      samp_plot <- NULL
      samp_plot_cov <- NULL
      }
    indiv_paths     <- coefs[, c("id", "lhs", "op", "rhs", "est.std", 
                                 "se", "z", "pvalue", "level")]
    indiv_paths$lhs <- recode.vars(indiv_paths$lhs, dat$lvarnames, dat$varnames)
    indiv_paths$rhs <- recode.vars(indiv_paths$rhs, dat$lvarnames, dat$varnames)
    indiv_paths     <- indiv_paths[order(indiv_paths$id, indiv_paths$level), ]
    colnames(indiv_paths) <- c("file", "lhs","op", "rhs", "beta", "se", 
                               "z", "pval", "level")
    # } # end "if no coefficients" commented out stl 11.20.17
    # combine fit information for summaryFit.csv
    
    fits        <- as.data.frame(do.call(rbind, store$fits))
    fits$file   <- rownames(fits)
    fits$status <- do.call(rbind, store$status)
    fits        <- fits[ ,c(12, 1:11, 13)]
    
    if (dat$subgroup){
      fits <- merge(fits, sub$sub_mem[ ,c(1,3)], by.x = "file", by.y = "names")  
      fits$modularity <- c(round(sub$modularity, digits = 4), 
                           rep("", (nrow(fits) - 1)))
      indiv_paths <- merge(indiv_paths, sub$sub_mem[ ,c(1,3)], 
                           by.x = "file", by.y = "names")
    }
    
    if (!is.null(dat$out)){ #& length(coefs[,1]) > 0){ # commented out stl 11.20.17
      write.csv(indiv_paths, file.path(dat$out, "indivPathEstimates.csv"),
                row.names = FALSE)
      write.csv(sample_counts, file.path(dat$out,
                                         "summaryPathCountsMatrix.csv"),
                row.names = FALSE)
      write.csv(fits, file.path(dat$out, "summaryFit.csv"), row.names = FALSE)
      if (dat$subgroup)
      write.csv(sub$sim, file.path(dat$out, "similarityMatrix.csv"), row.names = FALSE)
    }
    
  } else {
    indiv_paths <- store$coefs[[1L]]
    indiv_paths$file <- "all"
    indiv_paths$lhs  <- recode.vars(indiv_paths$lhs, dat$lvarnames, dat$varnames)
    indiv_paths$rhs  <- recode.vars(indiv_paths$rhs, dat$lvarnames, dat$varnames)
    indiv_paths      <- indiv_paths[order(indiv_paths$file), ]
    indiv_paths      <- indiv_paths[ ,c("file", "lhs","op", "rhs", "est.std", 
                                        "se", "z", "pvalue")]
    colnames(indiv_paths) <- c("file", "lhs", "op", "rhs", "beta", "se", "z", "pval")
    
    fits          <- store$fits[[1L]]
    file          <- c("all")
    names(file)   <- "file"
    status        <- store$status[[1L]]
    names(status) <- "status"
    fits          <- c(file, fits, status)
    fits <- t(fits)
    
    if (!is.null(dat$out)){
      write.csv(indiv_paths, file.path(dat$out, "allPathEstimates.csv"), 
                row.names = FALSE)
      write.csv(fits, file.path(dat$out, "summaryFit.csv"), row.names = FALSE)
    }
    
    sample_counts <- NULL
    sample_counts_corr <- NULL
    samp_plot_cov <- NULL
    samp_plot     <- NULL
  }
  
  dx <- list()
  if(diagnos){
    dx[[1]]<- dat
    dx[[2]] <- grp
    dx[[3]] <- store
    names(dx) <- c("dat", "grp", "store")
    }

    
  res <- list(fit           = fits,
              param_est     = indiv_paths,
              samp_plot     = samp_plot,
              samp_plot_cov = samp_plot_cov,
              sub_plots     = sub_plots,
              sub_plots_cov  = sub_plots_cov, 
              sample_counts = sample_counts,
              sample_counts_cov =    sample_counts_corr,
              sub_counts    = sub_counts,
              sub_counts_cov = sub_counts_cov,
              dx)  
  return(res)
  
}
