#' solution.tree
#' @name solution.tree
#' @title Solution trees for multiple solutions gimme.
#' @description This function allows for the exploration of divergences in multiple
#' solutions gimme for both the group and individuals.
#' @usage 
#' solution.tree(x,
#'               level     =  c("group", "individual"),
#'               cols      =  NULL,
#'               ids       =  "all",
#'               plot.tree =  FALSE)
#' @param x A fitted gimme object.
#' @param level A character vector indicating what levels of the solution tree 
#' you would like returned.  Options are "group", "individual", or c("group", "individual").
#' Defaults to c("group", "individual").
#' @param cols A character vector indicating additional information to include in tree plot.
#' Options include "stage", "pruned", "rmsea", "nnfi", "cfi","srmr", "grp_sol",
#' "bic", "aic", "modularity." Defaults to NULL.  
#' @param ids A character vector indicating the names of subjects to print.  Defaults to "all."
#' @param plot.tree Logical.  If TRUE, plot of tree is produced.  Defaults to FALSE.
#' @export solution.tree
solution.tree <- function(x, level = c("group", "individual"), cols = NULL, ids = "all", plot.tree = FALSE){
  
  tree <- batch.create.tree(
    x$grp_hist, 
    x$ind_hist, 
    x$ind_fit, 
    x$dat$subgroup, 
    names(x$dat$ts_list),
    x$sub)
  
  if(is.null(cols)){
    cols <- c("stage")
  } else {
    cols <- c("stage", cols)
  }
  
  if(level == "group"){
    
    grp_tree <- tree[[1]]
    
    if(plot.tree){
      final.plot <- plot(grp_tree)
      return(final.plot)
    } else {
      do.call(print, c(grp_tree, cols))
    }
    
  } else {

    ind_tree <-  tree[[2]]
    
    
    if(length(ids == 1) && ids == "all"){
      
      to_print <- seq(1:length(x$dat$ts_list))
      
    } else {
      
      to_print <- which(names(x$dat$ts_list) %in% ids)
      
    }
    
    if(length(ids == 1) && ids != "all"){
      
      if(plot.tree){
        plot(ind_tree[[to_print]])
      } else {
        #do.call(print, c(ind_tree, cols))
        print(c(ind_tree[[to_print]],cols), pruneMethod=NULL)
      }
      
    } else {
      
      for(i in to_print){
        cat(paste0(names(x$dat$ts_list)[i]),"\n\n")
        if(plot.tree){
          plot(ind_tree[[i]])
        } else {
          do.call(print, c(ind_tree[[i]], cols))
        }
      }
      
    }
    
    
  }
    
}
