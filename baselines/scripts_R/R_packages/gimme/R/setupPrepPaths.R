#' Allows user to open and close certain paths.
#' @param paths \code{lavaan}-style syntax containing paths with which
#' to begin model estimation (optional). That is, Y~X indicates that Y
#' is regressed on X, or X predicts Y. If no header is used,
#' then variables should be referred to with V followed (with no separation)
#' by the column number. If a header is used, variables should be referred to using 
#' variable names. To reference lag variables, "lag" should be added to the 
#' end of the variable name with no separation. Defaults to NULL.
#' @keywords internal
setupPrepPaths  <- function(paths, varLabels, ctrlOpts){
  
  # Satisfy CRAN checks
  varnames = NULL
  lvarnames = NULL
  
  table   <- lavaan::lavParTable(paths)
  
  # only include paths in the syntax which are specified free by the user
  # allows the possibility for the user to fix certain paths to zero
  tableFree    <- table[table$op == "~" & table$free != 0, ]
  
  dvsFree      <- tableFree$lhs
  ivsFree      <- tableFree$rhs
  
  # check if any exogenous variables have been incorrectly specified
  # for free paths
  if(!is.null(varLabels$uexo)){
    for (exog in varLabels$uexo){
      if (exog %in% dvsFree){
        stop(paste0('gimme ERROR: an exogenous variable was treated as endogenous in 
                    specified paths.  Please remove variable from exogenous list or 
                    correct path specification'))
      }
    }
  }
  
  if (nrow(tableFree) != 0){
    
    vsFree <- paste0(dvsFree, "~", ivsFree)
    
  } else {
    
    vsFree <- NULL
    
  }
  
  
  # table up the paths which are fixed to a certain value by the user
  tableFixed   <- table[table$op == "~" & table$free == 0,]
  
  if (nrow(tableFixed) > 0){
    
    dvsFixed     <- tableFixed$lhs
    
    # check if any exogenous variables have been incorrectly specified
    # for fixed paths
    if(!is.null(varLabels$uexo)){
      for (exog in varLabels$uexo){
        if (exog %in% dvsFixed){
          stop(paste0('gimme ERROR: an exogenous variable was treated as endogenous in 
                      specified paths.  Please remove variable from exogenous list or 
                      correct path specification'))
        }
      }
    }
    
    ivsFixed     <- recode.vars(tableFixed$rhs, varnames, lvarnames)
    vsFixed      <- paste0(dvsFixed, "~", ivsFixed)
    
  } else {
        
    vsFixed <- NULL
      
  }
  
  list = list(
    paths  = vsFree,
    remove = vsFixed
  )
  
  return(list)
}


