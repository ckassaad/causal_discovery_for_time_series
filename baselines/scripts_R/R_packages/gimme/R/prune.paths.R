#' Prunes paths. Ties together lowest.z and return.zs functions. 
#' @param base_syntax A character vector containing syntax that never changes.
#' @param fixed_syntax A character vector containing syntax that does not change
#' in a given stage of pruning.
#' @param add_syntax A character vector containing the syntax that is allowed
#' to change in a given stage of pruning.
#' @param data_list A list of datasets to be used in a given stage of the 
#' search. Varies based on group, subgroup, or individual-level stage.
#' @param n_paths The number of paths that are eligible for pruning. Equal
#' to the number of paths in add_syntax.
#' @param n_subj The number of subjects in a given stage of the search. If
#' in the group stage, n_subj equals the number of subjects. If in the subgroup
#' stage, n_subj equals the number of individuals in a given subgroup. At the 
#' individual stage, n_subj = 1.
#' @param prop_cutoff The proportion of individuals for whom a path must be
#' nonsignificant in order for it to be dropped from the models. NULL if used 
#' at the individual-level.
#' @param elig_paths A character vector containing eligible paths that
#' gimme is allowed to drop from the model at a given stage.
#' @param subgroup_stage Logical. Only present in order to instruct gimme
#' what message to print to console using writeLines.
#' @return Returns updated values of n_paths and add_syntax.
#' @keywords internal 
prune.paths <- function(base_syntax, 
                        fixed_syntax,
                        add_syntax,
                        data_list, 
                        n_paths, 
                        n_subj, 
                        prop_cutoff, 
                        elig_paths, 
                        subgroup_stage = FALSE){
  
  #-----------------------------------------------#
  # Determine the stage                           #
  #-----------------------------------------------#
  if(subgroup_stage){
    
    stage <- "subgroup prune" 
    
  } else {
    
    stage <- ifelse(n_subj == 1, "individual prune", "group prune")
    
  }
  #-----------------------------------------------#
  
  prune <- TRUE
  
  while(prune){      
    
    z_list <- list()
    
    for (k in 1:n_subj){
      
      if (!is.null(prop_cutoff)){
        
        if (is.list(fixed_syntax)){
          
          fit <- fit.model(
            syntax = c(base_syntax,fixed_syntax[[k]], add_syntax),
            data_file = data_list[[k]]
          )
          
        } else {
          
          fit <- fit.model(
            syntax = c(base_syntax,fixed_syntax, add_syntax),
            data_file = data_list[[k]]
          )
          
        }
        
        
        if (subgroup_stage){
          
          writeLines(paste0("subgroup-level pruning, subject ", k))
          
        } else {
          
          if(stage == "individual prune"){
            
            writeLines(paste0("individual-level pruning, subject ", k, " (",names(data_list)[k],")"))
            
          } else {
            
            writeLines(paste0("group-level pruning, subject ", k, " (",names(data_list)[k],")"))
            
          }
          
        }
        
      } else{
        
        fit <- fit.model(
          syntax = c(base_syntax, fixed_syntax, add_syntax),
          data_file = data_list
        ) 
        
      }
      
      z_list[[k]] <- return.zs(fit)
      
    }
    
    if(!all(is.na(z_list))){
      drop_param <- lowest.z(z_list,
                             elig_paths  = elig_paths,
                             prop_cutoff = prop_cutoff, 
                             n_subj      = n_subj)
    } else {
      drop_param <- NA
    }
    
    if (!is.na(drop_param)){
      n_paths    <- n_paths - 1
      add_syntax <- add_syntax[!add_syntax %in% drop_param] 
      
    } else {
      
      prune <- FALSE
      
    }
  }
  
  #---------------------------------------------------------#
  #---------------------------------------------------------#
  # syntax may be a list
  #---------------------------------------------------------#
  #---------------------------------------------------------#
  # if (is.list(fixed_syntax)){
  #   syntax <- lapply(seq_along(fixed_syntax), function(k){
  #     c(base_syntax, fixed_syntax[[k]], add_syntax)
  #   })
  # } else {
  #   syntax <-  c(base_syntax,fixed_syntax, add_syntax)
  # }
  # 
  # rmv_syntax <- setdiff(orig_add_syntax, add_syntax)
  # rmv_syntax <- if (length(rmv_syntax)==0) { rmv_syntax <- NA }
  
  res <- list(
    add_syntax = add_syntax,
    n_paths    = n_paths
  )
  
  return(res)
}