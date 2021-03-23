#' Identifies highest MI from list of MIs.
#' @param mi_list A list of MIs across individuals
#' @param indices A list of fit indices. Only relevant at the individual-level.
#' @param elig_paths A character vector containing eligible paths that
#' gimme is allowed to add to a model (e.g., no nonsense paths).
#' @param prop_cutoff The proportion of individuals for whom a path must be
#' significant in order for it to be added to the models. NULL if used 
#' at the individual-level.
#' @param n_subj The number of subjects in a given stage of the search. If
#' in the group stage, n_subj equals the number of subjects. If in the subgroup
#' stage, n_subj equals the number of individuals in a given subgroup. At the 
#' individual stage, n_subj = 1.
#' @param chisq_cutoff Cutoff used in order for MI to be considered significant.
#' Value varies depending on stage of search (e.g., group, subgroup, 
#' individual).
#' @return Returns name of parameter associated with highest MI. If no MI meets 
#' the criteria, returns NA.
#' @keywords internal 
highest.mi <- function(mi_list, 
                       indices,
                       elig_paths, 
                       prop_cutoff, 
                       n_subj,
                       chisq_cutoff,
                       allow.mult,
                       ms_tol,
                       hybrid,
                       dir_prop_cutoff){
  
  mi  = NULL # appease CRAN check
  sig = NULL # appease CRAN check
  param  = NULL # appease CRAN check
  
  mi_list_na       <- mi_list[!is.na(mi_list)]    # retain list form for later use
  n_converge    <- length(mi_list_na)
  mi_list       <- do.call("rbind", mi_list_na)
  
  mi_list$param <- paste0(mi_list$lhs, mi_list$op, mi_list$rhs)
  
  mi_list <- subset(mi_list, param %in% elig_paths, 
                    select = c("param", "mi", "epc"))
  
  mi_list$sig <- ifelse(mi_list$mi >= chisq_cutoff, 1, 0)
  
  mi_list <- transform(mi_list, 
                       sum = ave(mi, param, FUN = sum),
                       count = ave(sig, param, FUN = sum),
                       mean = ave(mi, param, FUN = mean))
  
  mi_list   <- subset(mi_list, !duplicated(param))
  mi_list   <- mi_list[order(-mi_list$count, -mi_list$sum), ]
  
  # we need to look at the means rather than the sum
  mi_list_ms <- mi_list[order(-mi_list$count, -mi_list$mean), ]
  
  #------------------------------------------------------#
  # Group search ongoing...
  #------------------------------------------------------#
  if (!is.null(prop_cutoff)){
    if(allow.mult){
      # if there are good solutions
      if (mi_list$count[1] > (prop_cutoff*n_converge)){
        
        red_mi     <- mi_list_ms[mi_list_ms$mean >= (mi_list_ms$mean[1]-ms_tol) & mi_list_ms$count == mi_list$count[1], , drop = FALSE]
        add_param  <- unlist(lapply(seq_along(1:nrow(red_mi)), function(i){
          red_mi$param[i]
        }))
        
        
      } else {
        
        add_param <- NA
        
      }
      
    } else {
      whichone <- 1
      go <- 0
      while(go == 0) { # for directionality search implemented after convo w Peter, Sy-Miin, & Jonathon
        
        mi_list_ms <- mi_list[order(-mi_list$count, -mi_list$mean), ]
        red_mi     <- mi_list_ms[whichone,]
        
        # if the opposite direction is significant for > proportion, go indiv by indiv and ensure its preferred in one direction > groupcutoff in one direction
        red_lhs <- strsplit(mi_list_ms$param[whichone], "~")[[1]][1]
        red_rhs <- strsplit(mi_list_ms$param[whichone], "~")[[1]][2]
        opposite <- paste0(red_rhs, "~", red_lhs)
        count_red <- 0
        count_opp <- 0 
        if (length(which(mi_list_ms$param == opposite))>0 && dir_prop_cutoff>0 && hybrid==FALSE) {
          if (!grepl("lag", mi_list_ms$param[whichone]) && mi_list[which(mi_list_ms$param == opposite), 6] >= (prop_cutoff*n_converge))
          {
            for (p in 1:length(mi_list_na)){
              grab_red_mi  <- mi_list_na[[p]][which(mi_list_na[[p]]$lhs == red_lhs),]
              red_mi_value <- grab_red_mi[which(grab_red_mi$rhs == red_rhs),]$mi
              grab_opp_mi  <- mi_list_na[[p]][which(mi_list_na[[p]]$lhs == red_rhs),]
              opp_mi_value <- grab_opp_mi[which(grab_opp_mi$rhs == red_lhs),]$mi
              ifelse(red_mi_value>opp_mi_value, (count_red = count_red+1), (count_opp = count_opp+1))
            }
          }
          if (((count_opp-count_red)/length(mi_list_na)) > dir_prop_cutoff){ 
            add_param<- opposite
            go <- 1} else if (((count_red-count_opp)/length(mi_list_na)) >= dir_prop_cutoff){
              add_param<- red_mi$param
              go <- 1} else{
                go <- 0 # directionality could not be determined for the group level, try next MI
                whichone <- whichone + 1}
        } else {
          go <- 1
          add_param <- ifelse(
            mi_list$count[whichone] > (prop_cutoff*n_converge),
            mi_list$param[whichone], 
            NA)
        }
        
      }
    }
    
    if (n_converge <= (n_subj/2)) { 
      
      add_param <- NA
    }
    #------------------------------------------------------#
    # Individual search ongoing...
    #------------------------------------------------------#
  } else {
    
    
    if(allow.mult){
      
      # we need to look at the means rather than the sum
      mi_list_ms <- mi_list[order(-mi_list$count, -mi_list$mean), ]
      red_mi     <- mi_list_ms[mi_list_ms$mean >= (mi_list_ms$mean[1]-ms_tol), , drop = FALSE]
      add_param  <- unlist(lapply(seq_along(1:nrow(red_mi)), function(i){
        red_mi$param[i]
      }))
      
    } else {
      
      add_param <- mi_list_ms$param[1L]
      
    }
    
    
    if (count.excellent(indices) >= 2) {
      
      add_param <- NA
      
    }
    
  }
  
  return(list(add_param=add_param,mi_list=mi_list))
}