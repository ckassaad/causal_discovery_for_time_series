#' Create structure of group search solutions.
#' @keywords internal   
subgroupStage <- function(dat,
                          grp,
                          confirm_subgroup,
                          elig_paths,
                          sub_feature,
                          sub_method,
                          ms_tol,
                          ms_allow,
                          sub_sim_thresh, 
                          hybrid,
                          dir_prop_cutoff){
# Satisfy CRAN checks
  sub_membership = NULL

    sub <- determine.subgroups(base_syntax  = c(dat$syntax, grp$group_paths),
                               data_list    = dat$ts_list,
                               n_subj       = dat$n_subj,
                               chisq_cutoff = dat$chisq_cutoff_mi_epc,
                               file_order   = dat$file_order,
                               elig_paths   = elig_paths,
                               confirm_subgroup = confirm_subgroup,
                               out_path     = dat$out, 
                               sub_feature  = sub_feature,
                               sub_method   = sub_method, 
                               sub_sim_thresh = sub_sim_thresh,
                               hybrid       = hybrid,
                               dir_prop_cutoff = dir_prop_cutoff)

  # begin subgroup-level search for paths ------------------------------------ #

  sub_spec <- vector("list", sub$n_subgroups)

  for (s in 1:sub$n_subgroups){

    sub_s <- list(sub_paths     = character(),
                  n_sub_paths   = 0,
                  sub_s_subjids = subset(sub$sub_mem,
                                         sub_membership == s)[ ,"names"],
                  n_sub_subj    = sum(sub$sub_mem$sub_membership == s,
                                      na.rm = TRUE),
                  sub_membership    = s)

    if (sub_s$n_sub_subj > 1){
      s4 <- search.paths(base_syntax  = dat$syntax,
                         fixed_syntax = grp$group_paths,
                         add_syntax   = character(),
                         n_paths      = 0,
                         data_list    = dat$ts_list[sub_s$sub_s_subjids],
                         elig_paths   = elig_paths,
                         prop_cutoff  = dat$sub_cutoff,
                         n_subj       = sub_s$n_sub_subj,
                         chisq_cutoff = qchisq(1-.05/sub_s$n_sub_subj, 1),
                         subgroup_stage = TRUE,
                         ms_tol         = ms_tol,
                         ms_allow       = FALSE)
      #sub_s[c("n_sub_paths", "sub_paths")] <- s4
      
      sub_s$sub_paths   <- s4[[1]][[1]]$add_syntax
      sub_s$n_sub_paths <- s4[[1]][[1]]$n_paths
      
    }
    sub_spec[[s]] <- sub_s
  }
  # end subgroup-level search for paths -------------------------------------- #

  # begin subgroup-level pruning --------------------------------------------- #
  for (s in 1:sub$n_subgroups){
    prune <- sub_spec[[s]]$n_sub_paths != 0 & sub_spec[[s]]$n_sub_subj != 1
    if(prune){
      s5 <- prune.paths(base_syntax  = dat$syntax,
                        fixed_syntax = grp$group_paths,
                        add_syntax   = sub_spec[[s]]$sub_paths,
                        data_list    = dat$ts_list[sub_spec[[s]]$sub_s_subjids],
                        n_paths      = sub_spec[[s]]$n_sub_paths,
                        n_subj       = sub_spec[[s]]$n_sub_subj,
                        prop_cutoff  = dat$sub_cutoff,
                        elig_paths   = sub_spec[[s]]$sub_paths,
                        subgroup_stage = TRUE)
      
      #sub_spec[[s]][c("n_sub_paths", "sub_paths")] <- s5
      sub_spec[[s]]$n_sub_paths <- s5$n_paths
      sub_spec[[s]]$sub_paths   <- s5$add_syntax
    }
  }

  # begin second-round group-level pruning ----------------------------------- #
  prune <- any(lapply(sub_spec, FUN = function(x) x$n_sub_paths != 0) == TRUE)

  sub_spec_comb <- do.call(rbind, sub_spec)
  ind           <- merge(sub$sub_mem, sub_spec_comb, "sub_membership", all.x = TRUE)
  ind           <- ind[order(ind$index),]
  ind$sub_paths[is.na(ind$sub_paths)] <- ""
  temp_count    <- grp$n_group_paths

  if (prune){
    s6 <- prune.paths(base_syntax  = dat$syntax,
                      fixed_syntax = ind$sub_paths,
                      add_syntax   = grp$group_paths,
                      data_list    = dat$ts_list,
                      n_paths      = grp$n_group_paths,
                      n_subj       = dat$n_subj,
                      prop_cutoff  = dat$group_cutoff,
                      elig_paths   = grp$group_paths,
                      subgroup_stage = FALSE)

    #grp[c("n_group_paths", "group_paths")] <- s6
    grp$n_group_paths <- s6$n_paths
    grp$group_paths   <- s6$add_syntax
  }

  if (temp_count != grp$n_group_paths){
    temp_sub_spec <- sub_spec
    for (s in 1:sub$n_subgroups){
      if (sub_spec[[s]]$n_sub_subj > 1){
        s7 <- search.paths(base_syntax  = dat$syntax,
                           fixed_syntax = grp$group_paths,
                           add_syntax   = sub_spec[[s]]$sub_paths,
                           n_paths      = sub_spec[[s]]$n_sub_paths,
                           data_list    =
                             dat$ts_list[sub_spec[[s]]$sub_s_subjids],
                           elig_paths   = elig_paths,
                           prop_cutoff  = dat$sub_cutoff,
                           n_subj       = sub_spec[[s]]$n_sub_subj,
                           chisq_cutoff =
                             qchisq(1-.05/sub_spec[[s]]$n_sub_subj, 1),
                           subgroup_stage = FALSE)
        #sub_spec[[s]][c("n_sub_paths", "sub_paths")] <- s7
        sub_spec[[s]]$sub_paths   <- s7[[1]][[1]]$add_syntax
        sub_spec[[s]]$n_sub_paths <- s7[[1]][[1]]$n_paths
      }
    }

    if (!identical(temp_sub_spec, sub_spec)){
      
      for (s in 1:sub$n_subgroups){
        prune <- temp_sub_spec[[s]]$n_sub_paths != sub_spec[[s]]$n_sub_paths
        if(prune){
          s8 <- prune.paths(base_syntax  = dat$syntax,
                            fixed_syntax = grp$group_paths,
                            add_syntax   = sub_spec[[s]]$sub_paths,
                            data_list    =
                              dat$ts_list[sub_spec[[s]]$sub_s_subjids],
                            n_paths      = sub_spec[[s]]$n_sub_paths,
                            n_subj       = sub_spec[[s]]$n_sub_subj,
                            prop_cutoff  = dat$sub_cutoff,
                            elig_paths   = sub_spec[[s]]$sub_paths,
                            subgroup_stage = FALSE)
          #sub_spec[[s]][c("n_sub_paths", "sub_paths")] <- s8
           sub_spec[[s]]$sub_paths   <- s8$add_syntax
           sub_spec[[s]]$n_sub_paths <- s8$n_paths
        }
      }
      
    }
    
  }

  sub_spec_comb <- do.call(rbind, sub_spec)
  ind           <- merge(sub$sub_mem, sub_spec_comb, "sub_membership", all.x = TRUE)
  ind$sub_paths[is.na(ind$sub_paths)] <- ""
  ind           <- ind[order(ind$index),]


  return(list(
    sub = sub,
    sub_spec = sub_spec,
    ind = ind,
    grp = grp
  ))
}