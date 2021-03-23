#' Searches for paths. Ties together highest.mi and return.mis functions.
#' @param base_syntax A character vector containing syntax that never changes.
#' @param fixed_syntax A character vector containing syntax that does not change
#' in a given stage of searching.
#' @param add_syntax A character vector containing the syntax that is allowed
#' to change in a given stage of searching.
#' @param data_list A list of datasets to be used in a given stage of the 
#' search. Varies based on group, subgroup, or individual-level stage.
#' @param n_paths The number of paths present in a given stage of searching.
#' Equal to the number of paths in add_syntax.
#' @param n_subj The number of subjects in a given stage of the search. If
#' in the group stage, n_subj equals the number of subjects. If in the subgroup
#' stage, n_subj equals the number of individuals in a given subgroup. At the 
#' individual stage, n_subj = 1.
#' @param prop_cutoff The proportion of individuals for whom a path must be
#' nonsignificant in order for it to be dropped from the models. NULL if used 
#' at the individual-level.
#' @param elig_paths A character vector containing eligible paths that
#' gimme is allowed to add to the model at a given stage.
#' @param chisq_cutoff Cutoff used in order for MI to be considered significant.
#' Value varies depending on stage of search (e.g., group, subgroup, 
#' individual).
#' @param subgroup_stage Logical. Only present in order to instruct gimme
#' what message to print to console using writeLines.
#' @return Returns updated values of n_paths and add_syntax.
#' @keywords internal 
search.paths.ms <- function(obj, 
                            data_list,
                            base_syntax, 
                            fixed_syntax,
                            elig_paths, 
                            prop_cutoff, 
                            n_subj, 
                            chisq_cutoff,
                            subgroup_stage,
                            ms_allow,
                            ms_tol,
                            hybrid,
                            dir_prop_cutoff){
  
  
  #-----------------------------------------------#
  # Determine the stage                           #
  #-----------------------------------------------#
  if(subgroup_stage){
    
    stage <- "subgroup search" 
    
  } else {
    
    stage <- ifelse(n_subj == 1, "individual search", "group search")
    
  }
  #-----------------------------------------------#
  
  
  new.obj <- list()
  
  for(j in 1:length(obj)){
    
    if(!obj[[j]]$final.sol){
      
      search    <- TRUE
      
      while(search){ # begin search
        
        mi_list <- list() 
        
        indices <- NULL
    
          # kmg 04.30.2019 remove for loop for individuals; use lapply 
        
          if (!is.null(prop_cutoff)){
            
            if(!subgroup_stage){
              writeLines(paste0("group-level search"))
            } else {
              writeLines(paste0("subgroup-level search"))
            }
           
            fit <- lapply(seq_along(data_list), function(i){fit.model(
              syntax= c(base_syntax, fixed_syntax, obj[[1]]$add_syntax),
              data_file = data_list[[i]])
            })

            for (k in 1:n_subj)
            mi_list[[k]] <- return.mis(fit[[k]], elig_paths)
          } else {
            
            for (k in 1:n_subj){
            
            # individual level search
            fit <- fit.model(
              syntax = c(base_syntax, fixed_syntax, obj[[j]]$add_syntax),
              data_file = data_list
            )
            
            #------------------------------------------------------#
            # Check to see if model converged.
            #------------------------------------------------------#
            
            if (!"error" %in% class(fit)){
              # stl 2018/08/16 separated convergence check from error check
              # can't inspect convergence of an error object
              if (lavaan::lavInspect(fit, "converged")){ 
                indices    <- fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", 
                                                 "srmr", "nnfi", "cfi"))
              } else indices <- NULL
            } else indices <- NULL
            mi_list[[k]] <- return.mis(fit, elig_paths)
            
          }
           }
        
        
        
        #------------------------------------------------------#
        # Add the parameters with the largest MI
        #------------------------------------------------------#
        if (!all(is.na(mi_list))){
          if (ms_allow | hybrid)
          add_p     <- highest.mi(mi_list      = mi_list,
                                  indices      = indices,
                                  elig_paths   = elig_paths,
                                  prop_cutoff  = prop_cutoff, 
                                  n_subj       = n_subj,
                                  chisq_cutoff = chisq_cutoff,
                                  allow.mult   = TRUE,
                                  ms_tol       = ms_tol,
                                  hybrid       = hybrid, 
                                  dir_prop_cutoff = dir_prop_cutoff)
          if (!ms_allow)
            add_p     <- highest.mi(mi_list      = mi_list,
                                    indices      = indices,
                                    elig_paths   = elig_paths,
                                    prop_cutoff  = prop_cutoff, 
                                    n_subj       = n_subj,
                                    chisq_cutoff = chisq_cutoff,
                                    allow.mult   = FALSE,
                                    ms_tol       = ms_tol,
                                    hybrid       = hybrid, 
                                    dir_prop_cutoff = dir_prop_cutoff)
          
          add_param <- add_p$add_param
          mi_info   <- add_p$mi_list
          
        } else {
          
          add_param            <- NA
          mi_info              <- NA
           
        }
        #------------------------------------------------------#
       
        
        #------------------------------------------------------#
        # If there are no paths to add.
        #------------------------------------------------------#
        if(all(is.na(add_param))){
          
          search               <- FALSE
          obj[[j]]$final.sol   <- TRUE
          
          res      <- list()
          res[[1]] <- list(
            add_syntax     = obj[[j]]$add_syntax,
            n_paths        = obj[[j]]$n_paths,
            final.sol      = obj[[j]]$final.sol
          )
        
        #------------------------------------------------------#
        # If there is only 1 path to add, 
        #  still searching...
        #------------------------------------------------------#
        } else if (length(add_param) == 1 | !ms_allow){
          
          
          obj[[j]]$n_paths     <- obj[[j]]$n_paths + 1
          obj[[j]]$add_syntax  <- append(obj[[j]]$add_syntax, add_param[1])
    
          
        #------------------------------------------------------#
        #  If there is > 1 path to add, stop searching.
        #------------------------------------------------------#  
        } else if (length(add_param) > 1 & ms_allow){
          
          writeLines(paste0("multiple solutions found..."))
          
          search <- FALSE # end the search
            
          n_paths    <- replicate(length(add_param), list(obj[[j]]$n_paths + 1))
          
          add_syntax <- replicate(length(add_param), obj[[j]]$add_syntax, simplify = FALSE)
          
          add_syntax <- lapply(seq_along(add_syntax), function(g) {
            append(obj[[j]]$add_syntax, add_param[g])
          })
          
  
          res <- list()
          
          for(i in 1:length(add_syntax)){
            
            res[[i]] <- list(
              add_syntax     = add_syntax[[i]],
              n_paths        = n_paths[[i]], 
              final.sol      = obj[[j]]$final.sol
            )
            
            
          }
            
        } 
        
      } # end search
      
      new.obj <- append(new.obj, res)
    
    } else {
    
      new.obj <- append(new.obj, obj[j])
    
    }
    
  }

  return(new.obj)
  
}