#' Set up base syntax file.
#' @keywords internal
setupBaseSyntax  <- function(paths, varLabels, ctrlOpts){
  
    #-------------------------------------------------------------#
    # NULL MODEL CONTAINS
    #-------------------------------------------------------------#
      # Cov & Var among variables than cannot be predicted.
      # Var of variables than can be predicted.
      # Means of exogenous variables
      # Intercepts of endogenous variables
      # Nonsense paths (fixed to zero)
      # Any fixed paths
    #-------------------------------------------------------------#
  
    # Var among variables than can be predicted.
    var.endo <- paste0(varLabels$endo, "~~", varLabels$endo)
  
    # Intercepts of endogenous variables
    int.endo  <- paste0(varLabels$endo, "~1")
  
    # Cov & Var among variables than cannot be predicted.
    
  ### Removes "conceptually exogenous" variables that are not statistically exogenous
  ### so covariances won't be estimated for those
    
  if(!is.null(varLabels$exog)){
    exog_model <- varLabels$exog[!varLabels$exog%in%varLabels$exog_lag]
  } else {
    exog_model <- varLabels$exog
  }
    
    
    cov.exog <- outer(exog_model, exog_model, function(x, y) paste0(x, "~~", y))
    cov.exog <- cov.exog[lower.tri(cov.exog, diag = TRUE)]
    
    # Means of exogenous variables
    mean.exog <- paste0(exog_model, "~1")
    
    # Nonsense paths (fixed to zero)
    nons.reg <- c(t(outer(varLabels$exog, varLabels$endo, function(x, y) paste0(x, "~0*", y))))
   
    # Nonsense paths (not fixed to zero)
    nons.paths <- c(t(outer(varLabels$exog, varLabels$endo, function(x, y) paste0(x, "~", y))))
    
    # Any fixed paths
    fixed.paths <- paths
    
    if(ctrlOpts$ar) {
      ar.paths    <- paste0(
        setdiff(varLabels$orig, varLabels$exog_con), 
        "~", paste0(setdiff(varLabels$orig, varLabels$exog_con),"lag")
      )
      fixed.paths <- c(fixed.paths, ar.paths)
    }
    
    # All Possible Paths
    all.poss <- outer(varLabels$endo, c(varLabels$endo, varLabels$exog), function(x, y) paste0(x, "~", y))
    all.poss <- c(all.poss[lower.tri(all.poss, diag = FALSE)], all.poss[upper.tri(all.poss, diag = FALSE)])
    
    # All Possible Contemporaneous Correlations
    all.corr <- outer(varLabels$endo, varLabels$endo, function(x,y) paste0(x, "~~", y))
    all.corr <- c(all.corr[lower.tri(all.corr, diag = FALSE)], all.corr[upper.tri(all.corr, diag = FALSE)])
    # Both V1~~V2 and V2~~V1 are kept for now because we don't know which one lavaan produces
    
     # If multiplied vars are specified, remove prediction of first order effects by multiplied
     # variables
     nons.mult<-character(2*length(varLabels$mult))
     j<-1
     if(!is.null(varLabels$mult)){
       for(i in 1:length(varLabels$mult)){
          mult.names <- unlist(strsplit(varLabels$mult[i], "_by_", fixed = TRUE))
          nons<-apply(expand.grid(mult.names[1:length(mult.names)],
                                  varLabels$mult[i]), 1, paste, collapse = "~")
 
           nons.mult[j]<-nons[1]
          nons.mult[j+1]<-nons[2]
          j<-(j+2)
     }
     }
 
      all.poss<-setdiff(all.poss, nons.mult)
    
    base.syntax <- c(
      var.endo,
      int.endo,
      cov.exog,
      mean.exog,
      nons.reg,
      fixed.paths
    )
    
    return(list(
      syntax = base.syntax,
      fixed.paths = fixed.paths,
      candidate.paths = setdiff(all.poss, fixed.paths),
      candidate.corr = all.corr,
      nonsense.paths =  nons.paths 
    ))
  
   
}


