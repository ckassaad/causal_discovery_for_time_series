#' Identifies lowest z value from list of z values.
#' @param z_list A list of z values across individuals.
#' @param elig_paths A character vector containing eligible paths that
#' gimme is allowed to drop from the model at a given stage.
#' @param prop_cutoff The proportion of individuals for whom a path must be
#' nonsignificant in order for it to be dropped from the models. NULL if used 
#' at the individual-level.
#' @param n_subj The number of subjects in a given stage of the search. If
#' in the group stage, n_subj equals the number of subjects. If in the subgroup
#' stage, n_subj equals the number of individuals in a given subgroup. At the 
#' individual stage, n_subj = 1.
#' @return Returns name of parameter associated with lowest z. If no z meets 
#' the criteria, returns NA.
#' @keywords internal 
lowest.z <- function(z_list, 
                     elig_paths, 
                     prop_cutoff,
                     n_subj){
  
  param  = NULL # appease CRAN check
  z      = NULL # appease CRAN check
  sig    = NULL # appease CRAN check
  
  z_list     <- z_list[!is.na(z_list)]
  n_converge <- length(z_list)
  z_list     <- do.call("rbind", z_list)
  z_list$param <- paste0(z_list$lhs, z_list$op, z_list$rhs)
  z_list       <- subset(z_list, param %in% elig_paths,
                         select = c("param", "z"))
  z_list$sig   <- ifelse(abs(z_list$z) > abs(qnorm(.025/n_subj)), 1, 0)
  z_list <- transform(z_list,
                      sum   = ave(abs(z), param, FUN = sum),
                      count = ave(sig, param, FUN = sum))
  z_list <- subset(z_list, !duplicated(param))
  z_list <- z_list[order(z_list$count, z_list$sum), ]
  if (!is.null(prop_cutoff)){
    drop_param <- ifelse(z_list$count[1L] <= (prop_cutoff*n_converge),
                         z_list$param[1L], NA)
    if (n_converge <= (n_subj/2)) drop_param <- NA
  } else {
    drop_param <- ifelse(z_list$sig[1L] == 0, z_list$param[1L], NA)
  }
  return(drop_param)
}
