#' @name convolve
#' @aliases  convolve setupConvolve
#' @title Group iterative multiple model estimation.
#' @description This function estimates the basis vectors related to responses following 
#' a binary impulse and convolves that binary impulse vector. 
#' @usage
#' convolveFIR(ts_list = NULL, 
#'      varLabels = NULL, 
#'      conv_length = 16, 
#'      conv_interval = 1)
#'
#' @param ts_list a list of dataframes.
#' @param varLabels a list of variable sets. Contains varLabels$coln, all column names, varLabels$conv, 
#' the names of variables to convolve, and varLabels$exog, a list of exogenous variables (if any).
#' @param conv_length Expected response length in seconds. For functional MRI BOLD, 16 seconds (default) is typical
#' for the hemodynamic response function. 
#' @param conv_interval Interval between data acquisition. Currently must be a constant. For 
#' fMRI studies, this is the repetition time. Defaults to 1. 
#' @keywords setupConvolve
#' @export convolveFIR
convolveFIR <- setupConvolve <- function(ts_list = NULL, 
                          varLabels = NULL, 
                          conv_length = 16, 
                          conv_interval = 1){
  
  # Satisfy CRAN checks
 ts = NULL
  # We only convolve contemporaneous (lagged contemporaneous created afterwards). 
  to_convolve <- setdiff(varLabels$coln, c(varLabels$conv, varLabels$exog))
  
  ts_list <- lapply(ts_list, function(df){
    
    conv_use  <- df[,to_convolve, drop = FALSE]
    
    if(any(apply(conv_use, 2, function(x) any(is.na(x) | is.infinite(x))))){
      
      conv_use[]  <- apply(conv_use, 2, function(x) { imputeTS::na_kalman(stats::ts(x)) })
      
    }
    
    for (cv in varLabels$conv){
      
      stimuli   <- df[,cv, drop = TRUE]
      
      if(any(is.na(stimuli))){
        stop(
          "gimme ERROR: missing values in the binary impulse vector not allowed"
        )
      }
      
      convolved <- sFIR(data = conv_use, stimuli = stimuli, response_length = conv_length, interval = conv_interval)
      
      df[,cv]   <- convolved$conv_stim_onsets[1:nrow(df)]

    }
    
    df

  })
    
  return(ts_list)
  
}