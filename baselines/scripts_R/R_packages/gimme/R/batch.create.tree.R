#' Create tree structures for group search solutions.
#' @keywords internal 
batch.create.tree <- function(hist, ind_hist, ind_fit, subgroup, names.ts_list, sub){
  ######
  # hist <- x$grp_hist
  # ind_hist <- x$ind_hist
  # ind_fit <- x$ind_fit
  # subgroup <- x$dat$subgroup
  # names.ts_list <- names(x$dat$ts_list)
  # sub <- x$sub
  #####

  #-----------------------------------------------#
  # If subgroup == FALSE                          #
  #-----------------------------------------------#
  
  if(!subgroup){
    
    # ind_hist <- ind_hist_old; j <- 1
    final_grp_hist<- hist[[length(hist)]]
    final_grp_hist <- lapply(seq_along(final_grp_hist), function(i){
         final_grp_hist[[i]]$grp_sol <- i
         final_grp_hist[[i]]$fit     <- NA
         final_grp_hist[[i]]
    })

    trees_ind <- lapply(seq_along(1:length(ind_hist[[1]])), function(j){
      ind_hist_k <- lapply(seq_along(final_grp_hist), function(i){
        lapply(seq_along(ind_hist[[i]][[j]]), function(l){
          list(
            list(
              "subj" = names.ts_list[j],
              "add_syntax" =  final_grp_hist[[i]]$add_syntax,
              "group_sol" = ind_hist[[i]][[j]][[l]]$grp_sol,
              "ind_sol" = NA,
              "pruned" = final_grp_hist[[i]]$pruned,
              "stage"  = "group",
              "fit" = NA
            ),
            list(
              "subj" = names.ts_list[j],
              "add_syntax" =  c(final_grp_hist[[i]]$add_syntax, ind_hist[[i]][[j]][[l]]$add_syntax),
              "group_sol" = ind_hist[[i]][[j]][[l]]$grp_sol,
              "ind_sol" = paste0(ind_hist[[i]][[j]][[l]]$grp_sol, "_",
                  ind_hist[[i]][[j]][[l]]$ind_sol),
              "pruned" =  ind_hist[[i]][[j]][[l]]$pruned,
              "stage"  = "individual",
              "fit" = ind_fit[[i]][[j]][[l]]$fits
            )
          )
        })
      })
      tree_ind <- create.tree(ind_hist_k, subgroup, individual = TRUE)
      tree_ind
    })
    
    
    if(final_grp_hist[[1]]$n_paths == 0){
      trees_grp <- NA
    } else {
      grp_hist_i <- lapply(seq_along(final_grp_hist), function(i){
        list(
          "subj"= NA,
          "add_syntax" =  final_grp_hist[[i]]$add_syntax,
          "grp_sol" = i,
          "ind_sol" = NA,
          "pruned" = final_grp_hist[[i]]$pruned,
          "stage"  = "group",
          "fit" = NA
        )
      })
      trees_grp <- create.tree(grp_hist_i, subgroup, individual = FALSE)
    }
      
    
    
    
    # ind_all <- lapply(seq_along(ind_hist), function(i){
    #   lapply(seq_along(ind_hist[[i]]), function(j){
    #     lapply(seq_along(ind_hist[[i]][[j]]), function(l){
    #       list(
    #         list(
    #           "subj_id" =  j,
    #           "add_syntax" =  final_grp_hist[[i]]$add_syntax,
    #           "group_sol"  = ind_hist[[i]][[j]][[l]]$grp_sol,
    #           "ind_sol"    = NA,
    #           "pruned"     = final_grp_hist[[i]]$pruned,
    #           "stage"      = "group"
    #         ),
    #         list(
    #           "subj_id" =  j,
    #           "add_syntax" =  c(final_grp_hist[[i]]$add_syntax, ind_hist[[i]][[j]][[l]]$add_syntax),
    #           "group_sol" = ind_hist[[i]][[j]][[l]]$grp_sol,
    #           "ind_sol" = ind_hist[[i]][[j]][[l]]$ind_sol,
    #           "pruned" =  ind_hist[[i]][[j]][[l]]$pruned,
    #           "stage"  = "individual"
    #         )
    #       )
    #     })
    #   })
    # })

    return(list(trees_grp = trees_grp, trees_ind = trees_ind))
    
  #-----------------------------------------------#
  # If subgroup == TRUE                           #
  #-----------------------------------------------#  
  } else if(subgroup){
    
    #-----------------------------------------------#
    # Individuals Plotted Individually              #
    #-----------------------------------------------#  
    
    final_grp_hist <- hist[[length(hist)-1]]
    final_sub_hist <- hist[[length(hist)]]
      
    trees <- lapply(seq_along(1:length(ind_hist[[1]])), function(j){
      
      ind_hist_k <- lapply(seq_along(final_grp_hist), function(i){
        
        subgrp_m  <- sub[[i]]$sub_mem[sub[[i]]$sub_mem$names == names.ts_list[j],"sub_membership"]
        grp_i_idx <- unlist(lapply(final_sub_hist, "[", c("grp_sol"))) == i
        sub_m_idx <- unlist(lapply(final_sub_hist, "[", c("sub_sol"))) == subgrp_m
        grp_i_sub_m_idx <- which(grp_i_idx & sub_m_idx)
        
        lapply(seq_along(ind_hist[[i]][[j]]), function(l){
          list(
            list(
              "subj" = names.ts_list[j],
              "add_syntax" =  final_grp_hist[[i]]$add_syntax,
              "grp_sol" = ind_hist[[i]][[j]][[l]]$grp_sol,
              "sub_sol" = ind_hist[[i]][[j]][[l]]$sub_sol,
              "ind_sol" = paste0(ind_hist[[i]][[j]][[l]]$grp_sol, "_",ind_hist[[i]][[j]][[l]]$ind_sol),
              "pruned" = final_grp_hist[[i]]$pruned,
              "stage"  = "group",
              "modularity" = NA,
              "fit" = NA
            ),
            list(
              "subj" = names.ts_list[j],
              "add_syntax" =  c(final_grp_hist[[i]]$add_syntax, final_sub_hist[[grp_i_sub_m_idx]]$add_syntax), 
              "grp_sol" = ind_hist[[i]][[j]][[l]]$grp_sol,
              "sub_sol" = ind_hist[[i]][[j]][[l]]$sub_sol,
              "ind_sol" = paste0(ind_hist[[i]][[j]][[l]]$grp_sol, "_",ind_hist[[i]][[j]][[l]]$ind_sol),
              "pruned" =  final_sub_hist[[grp_i_sub_m_idx]]$pruned,
              "stage"  = "subgroup",
              "modularity" = NA,
              "fit" = NA
            ),
            list(
              "subj" = names.ts_list[j],
              "add_syntax" =  c(final_grp_hist[[i]]$add_syntax, 
                final_sub_hist[[grp_i_sub_m_idx]]$add_syntax,
                ind_hist[[i]][[j]][[l]]$add_syntax
              ),
              "group_sol" = ind_hist[[i]][[j]][[l]]$grp_sol,
              "sub_sol" = ind_hist[[i]][[j]][[l]]$sub_sol,
              "ind_sol" = paste0(ind_hist[[i]][[j]][[l]]$grp_sol, "_",ind_hist[[i]][[j]][[l]]$ind_sol),
              "pruned" =  c(final_grp_hist[[i]]$pruned, ind_hist[[i]][[j]][[l]]$pruned),
              "stage"  = "individual",
              "modularity" = NA,
              "fit" = ind_fit[[i]][[j]][[l]]$fits
            )
          )
        })
      })
      trees_ind <- create.tree(ind_hist_k, subgroup, individual = TRUE)
    })
    
    
    
    #-----------------------------------------------#
    # Individuals Plotted Together                  #
    #-----------------------------------------------#  
    
    # final_grp_hist <- hist[[length(hist)-1]]
    # final_sub_hist <- hist[[length(hist)]]
    #   
    # all_ind <- lapply(seq_along(1:length(ind_hist[[1]])), function(j){
    #   
    #   ind_hist_k <- lapply(seq_along(final_grp_hist), function(i){
    #     
    #     subgrp_m  <- sub[[i]]$sub_mem[sub[[i]]$sub_mem$names == names.ts_list[j],"sub_membership"]
    #     grp_i_idx <- unlist(lapply(final_sub_hist, "[", c("grp_sol"))) == i
    #     sub_m_idx <- unlist(lapply(final_sub_hist, "[", c("sub_sol"))) == subgrp_m
    #     grp_i_sub_m_idx <- which(grp_i_idx & sub_m_idx)
    #     
    #     lapply(seq_along(ind_hist[[i]][[j]]), function(l){
    #       list(
    #         list(
    #           "subj" = names.ts_list[j],
    #           "add_syntax" =  final_grp_hist[[i]]$add_syntax,
    #           "grp_sol" = ind_hist[[i]][[j]][[l]]$grp_sol,
    #           "sub_sol" = ind_hist[[i]][[j]][[l]]$sub_sol,
    #           "ind_sol" = paste0(ind_hist[[i]][[j]][[l]]$grp_sol, "_",ind_hist[[i]][[j]][[l]]$ind_sol),
    #           "pruned" = final_grp_hist[[i]]$pruned,
    #           "stage"  = "group",
    #           "modularity" = NA,
    #           "fit" = NA
    #         ),
    #         list(
    #           "subj" = names.ts_list[j],
    #           "add_syntax" =  c(final_grp_hist[[i]]$add_syntax, final_sub_hist[[grp_i_sub_m_idx]]$add_syntax), 
    #           "grp_sol" = ind_hist[[i]][[j]][[l]]$grp_sol,
    #           "sub_sol" = ind_hist[[i]][[j]][[l]]$sub_sol,
    #           "ind_sol" = paste0(ind_hist[[i]][[j]][[l]]$grp_sol, "_",ind_hist[[i]][[j]][[l]]$ind_sol),
    #           "pruned" =  final_sub_hist[[grp_i_sub_m_idx]]$pruned,
    #           "stage"  = "subgroup",
    #           "modularity" = NA,
    #           "fit" = NA
    #         ),
    #         list(
    #           "subj" = names.ts_list[j],
    #           "add_syntax" =  c(final_grp_hist[[i]]$add_syntax, 
    #             final_sub_hist[[grp_i_sub_m_idx]]$add_syntax,
    #             ind_hist[[i]][[j]][[l]]$add_syntax
    #           ),
    #           "group_sol" = ind_hist[[i]][[j]][[l]]$grp_sol,
    #           "sub_sol" = ind_hist[[i]][[j]][[l]]$sub_sol,
    #           "ind_sol" = paste0(ind_hist[[i]][[j]][[l]]$grp_sol, "_",ind_hist[[i]][[j]][[l]]$ind_sol),
    #           "pruned" =  c(final_grp_hist[[i]]$pruned, ind_hist[[i]][[j]][[l]]$pruned),
    #           "stage"  = "individual",
    #           "modularity" = NA,
    #           "fit" = ind_fit[[i]][[j]][[l]]$fits
    #         )
    #       )
    #     })
    #   })
    # })
    
    
    #-----------------------------------------------#
    # Group subgroup History                        #
    #-----------------------------------------------#  
    
    final_grp_hist <- hist[[length(hist)-1]]
    final_sub_hist <- hist[[length(hist)]]
  
    grpsub_hist <- lapply(seq_along(final_sub_hist), function(i){
      list(
        list(
          "add_syntax" = final_grp_hist[[final_sub_hist[[i]]$grp_sol]]$add_syntax,
          "grp_sol"    = final_sub_hist[[i]]$grp_sol,
          "sub_sol"    = NA,
          "pruned"     = final_grp_hist[[final_sub_hist[[i]]$grp_sol]]$pruned,
          "stage"      = "group",
          "modularity" = NA,
          "fit"        = NA
        ),
        list(
          "add_syntax" = c(final_grp_hist[[final_sub_hist[[i]]$grp_sol]]$add_syntax,
             final_sub_hist[[i]]$add_syntax),              
          "grp_sol"    = final_sub_hist[[i]]$grp_sol,
          "sub_sol"    = final_sub_hist[[i]]$sub_sol,
          "pruned"     = final_sub_hist[[i]]$pruned,
          "stage"      = "subgroup",
          "modularity" = sub[[final_sub_hist[[i]]$grp_sol]]$modularity,
          "fit" = NA
        )
      )
    })
 
    trees_grpsub <- create.tree(grpsub_hist, subgroup, individual = FALSE)
    
   return(list(trees_grp = trees_grpsub, trees_ind = trees))
    
  }
}