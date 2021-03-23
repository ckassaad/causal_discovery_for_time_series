#' Individual-level search. Used in gimmeSEM, aggSEM, indSEM.
#' @param dat A list containing information created in setup().
#' @param grp A list containing group-level information. NULL in aggSEM and
#' indSEM.
#' @param ind A list containing individual- and (potentially) subgroup-level
#' information.
#' @return Lists associated with coefficients, fit indices, etc.
#' @keywords internal 
indiv.search.ms <- function(dat, grp, ind, ms_tol, ms_allow, grp_num, hybrid){
  
  #-----------------------------------------------#
  # Prepare data depending on proc (agg vs gimme) #
  #-----------------------------------------------#
  if (dat$agg){
    
    ind   <- NULL
    n_ind <- 1
    data_all <- do.call(rbind, dat$ts_list)
    colnames(data_all) <- dat$varnames
    
  } else {
    
    ind$ind_paths   <-  vector("list", dat$n_subj)
    ind$n_ind_paths <- 0
    n_ind           <- dat$n_subj
    
  }
  
  
  ind_all         <- replicate(n_ind, list(), simplify = FALSE)
  ind_history_all <- replicate(n_ind, list(), simplify = FALSE)
  
  #-----------------------------------------------#
  # Run individual search                         #
  #-----------------------------------------------#
  for (k in 1:n_ind){ # k <- 1
    
    if (dat$agg){
      
      data_list <- data_all
      
    } else {
      
      data_list <- dat$ts_list[[k]]
      
      writeLines(paste0("individual-level search, subject ", k, " (", names(dat$ts_list)[k],")"))
      
    }
    
    #-----------------------------------------------#
    # Conduct the ind level search                  #
    #-----------------------------------------------#
    s1.ind <- search.paths(
      base_syntax  = dat$syntax, 
      fixed_syntax = c(grp$group_paths, ind$sub_paths[[k]]),
      add_syntax   = character(),
      n_paths      = 0,
      data_list    = data_list,
      elig_paths   = dat$candidate_paths,
      prop_cutoff  = NULL,
      n_subj       = 1,
      chisq_cutoff = qchisq(.99, 1),
      subgroup_stage = FALSE,
      ms_allow       = TRUE,
      ms_tol         = ms_tol,
      hybrid         = hybrid
    )
    
    s1.ind <- lapply(seq_along(s1.ind), function(i){
      lapply(seq_along(s1.ind[[i]]), function(j){
        s1.ind[[i]][[j]]$pruned     <- NA
        s1.ind[[i]][[j]]$stage      <- "individual"
        s1.ind[[i]][[j]]$grp_sol    <- grp_num
        s1.ind[[i]][[j]]$sub_sol    <- ind$sub_membership[k]
        s1.ind[[i]][[j]]$ind_sol    <- NA
        if(i == length(s1.ind)){ s1.ind[[i]][[j]]$ind_sol <- j }
        s1.ind[[i]][[j]]
      })
    })
  

    #-----------------------------------------------#
    # Conduct the first pruning                     #
    #-----------------------------------------------#
    p1.ind <- lapply(seq_along(s1.ind[[length(s1.ind)]]), function(i){
      
      if(ms_allow){
        writeLines(paste0("pruning solution ", i ," for subject ", k, "..."))
      }
      
      prune.paths(
        base_syntax    = dat$syntax,
        fixed_syntax   = c(grp$group_paths, ind$sub_paths[[k]]),
        add_syntax     = s1.ind[[length(s1.ind)]][[i]]$add_syntax,
        data_list      = data_list,
        n_paths        = s1.ind[[length(s1.ind)]][[i]]$n_paths,
        n_subj         = 1,
        prop_cutoff    = NULL,
        elig_paths     = s1.ind[[length(s1.ind)]][[i]]$add_syntax,
        subgroup_stage = FALSE
      )
      
    })
    
    #-----------------------------------------------#
    # Add the pruning step to the indiv. history    #
    #-----------------------------------------------#
    
    s1.ind[[length(s1.ind)]] <- lapply(seq_along(s1.ind[[length(s1.ind)]]), function(i){
      pruned <- setdiff(s1.ind[[length(s1.ind)]][[i]]$add_syntax,p1.ind[[i]]$add_syntax)
      if(length(pruned) != 0){ s1.ind[[length(s1.ind)]][[i]]$pruned <- pruned} 
      s1.ind[[length(s1.ind)]][[i]]$stage      <- "individual"
      s1.ind[[length(s1.ind)]][[i]]$grp_sol    <- grp_num
      s1.ind[[length(s1.ind)]][[i]]$sub_sol    <- ind$sub_membership[k]
      s1.ind[[length(s1.ind)]][[i]]$ind_sol    <- s1.ind[[length(s1.ind)]][[i]]$ind_sol 
      s1.ind[[length(s1.ind)]][[i]]
    })
    
    
    #-----------------------------------------------#
    # One final search phase?                       #
    #-----------------------------------------------# 
    
    s2.ind <- sapply(seq_along(p1.ind), function(i){
      
      if(ms_allow){
        writeLines(paste0("searching solution ", i ," for subject ", k, "..."))
      }
      
      search.paths(
        base_syntax  = dat$syntax, 
        fixed_syntax = c(grp$group_paths, ind$sub_paths[[k]]),
        add_syntax   = p1.ind[[i]]$add_syntax,
        n_paths      = p1.ind[[i]]$n_paths,
        data_list    = data_list,
        elig_paths   = dat$candidate_paths,
        prop_cutoff  = NULL,
        n_subj       = 1,
        chisq_cutoff = 0,
        subgroup_stage = FALSE,
        ms_allow       = FALSE, # do not allow multiple solutions on cleanup
        ms_tol         = ms_tol,
        hybrid         = hybrid
      )
      
    }, simplify = TRUE)
    
    #-----------------------------------------------#
    # Construct final individual search history     #
    #-----------------------------------------------# 
  
    s1.ind[[length(s1.ind)]] <- lapply(seq_along(s1.ind[[length(s1.ind)]]), function(i){
      # are there any new paths to add
      new_paths <- setdiff(s2.ind[[i]][[1]]$add_syntax,s1.ind[[length(s1.ind)]][[i]]$add_syntax)
      if(length(new_paths) > 0){ 
        s1.ind[[length(s1.ind)]][[i]]$add_syntax <- s2.ind[[i]][[1]]$add_syntax
      }
      # do we have to add any pruned path added back in?
      if(!any(is.na(s1.ind[[length(s1.ind)]][[i]]$pruned))){
        prunes_left <- setdiff(s1.ind[[length(s1.ind)]][[i]]$pruned, s2.ind[[i]][[1]]$add_syntax)
        if(length(prunes_left) > 0){
          s1.ind[[length(s1.ind)]][[i]]$pruned <- prunes_left
        } else {
          s1.ind[[length(s1.ind)]][[i]]$pruned <- NA
        }
      }
      s1.ind[[length(s1.ind)]][[i]]
    })
    

    
    
    ind_k <- replicate(length(s1.ind[[length(s1.ind)]]), ind, simplify = FALSE)
    
    ind_k <- lapply(seq_along(ind_k), function(r){
       ind_k[[r]]$ind_paths[[k]] <- s1.ind[[length(s1.ind)]][[r]]$add_syntax
       ind_k[[r]]$n_ind_paths[k] <- s1.ind[[length(s1.ind)]][[r]]$n_paths
       ind_k[[r]]
    })
    
    ind_all[[k]] <- ind_k
    ind_history_all[[k]] <- s1.ind[[length(s1.ind)]] # just added
    
  }
  
    
  #-----------------------------------------------#
  # Prepare output                                #
  #-----------------------------------------------#
      
  res <- lapply(seq_along(1:n_ind), function(k){
    lapply(seq_along(ind_all[[k]]), function(j){
      s10 <- get.params(dat, grp, ind_all[[k]][[j]], k)
      list(
        "subj"   = names(dat$ts_list)[k],
        "ind"    = ind_all[[k]][[j]],
        "status" = s10$status,
        "fits"   = s10$ind_fit,
        "coefs"  = s10$ind_coefs,
        "betas"  = s10$ind_betas,
        "psi"    = s10$ind_psi, 
        "psiunstd"  = s10$ind_psi_unstd, 
        "vcov"   = s10$ind_vcov,
        "plots"  = s10$ind_plot,
        "syntax" = s10$ind_syntax
      )
    })
  })
      
      

  return(list(res = res, ind_history = ind_history_all))
}