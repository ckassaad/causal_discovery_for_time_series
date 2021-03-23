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
search.paths <- function(base_syntax, 
                         fixed_syntax,
                         add_syntax,
                         n_paths, 
                         data_list, 
                         elig_paths, 
                         prop_cutoff, 
                         n_subj, 
                         chisq_cutoff,
                         subgroup_stage = FALSE,
                         ms_allow = FALSE,
                         ms_tol   = 1e-6,
                         hybrid = F,
                         dir_prop_cutoff = 0){
  
  ####################################
  # ind search
  ####################################
  # base_syntax  = dat$syntax
  # fixed_syntax = c(grp$group_paths, ind$sub_paths[[k]])
  # add_syntax   = character()
  # n_paths      = 0
  # data_list    = data_list
  # elig_paths   = dat$candidate_paths
  # prop_cutoff  = NULL
  # n_subj       = 1
  # chisq_cutoff = qchisq(.99, 1)
  # subgroup_stage = FALSE
  # ms_allow       = TRUE
  # ms_tol         = 1e-5
  ####################################
  
  ####################################
  # group search
  ####################################
    # base_syntax    = dat$syntax
    # fixed_syntax   = NULL
    # add_syntax     = grp$group_paths
    # n_paths        = grp$n_group_paths
    # data_list      = dat$ts_list
    # elig_paths     = dat$candidate_paths
    # prop_cutoff    = dat$group_cutoff
    # n_subj         = dat$n_subj
    # chisq_cutoff   = qchisq(1-.05/dat$n_subj, 1)
    # subgroup_stage = FALSE
    # ms_allow       = ms_allow
    # ms_tol         = ms_tol
  ####################################

  obj <- replicate(1, 
    list(
      add_syntax     = add_syntax,
      n_paths        = n_paths, 
      final.sol      = FALSE
    ), simplify = FALSE
  )
  
  cnt <- 1
  history <- list()
  
  while(!all(unlist(lapply(obj,"[[","final.sol")))){
    
    obj <- search.paths.ms(
      obj, 
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
      dir_prop_cutoff
    )
    
    history[[cnt]] <- obj
    cnt            <- cnt + 1
    
  }
  
  # remove final solution from obj
  history <- lapply(history, function(obj){
     lapply(obj, "[", c("add_syntax", "n_paths"))
  })
  
  return(history)
  
 
}

