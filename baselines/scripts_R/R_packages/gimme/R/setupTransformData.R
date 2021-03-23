#' Transform raw data as required.
#'
#' @param ts_list a list or directory
#' @param varLabels Variable labels. 
#' @param ctrlOpts List used in setup function.
#' @keywords internal
setupTransformData <- function(ts_list       = NULL, 
                               varLabels     = NULL,
                               ctrlOpts      = NULL ){
  
  
  #-------------------------------------------------------------#
  # Early data checks (need to happen before any processing)
  #-------------------------------------------------------------# 
  
  # check to make sure variables in exogenous argument exist in data
  if(!is.null(varLabels$exog)){
    
    for(exog in setdiff(varLabels$exog, c(varLabels$mult, varLabels$lagg))){
      if (!exog %in% colnames(ts_list[[1]])){
        stop(paste0('gimme ERROR: Exogenous variable name provided is not in data column names
                    Please fix.'))
      }
    }
    
  }
  
  if (ncol(ts_list[[1]]) == 1) {
    
    stop(paste0("gimme ERROR: only one column of data read in. ",
                "Check if sep argument properly specified."))
    
  }
  
  
  #-------------------------------------------------------------#
  # (1) convolve variables at t.
  #-------------------------------------------------------------#
  # If neccessary, convolve  variables at t.
  #  * Missing values in variables used to generate convolution
  #    parameters imputed using imputeTS::na.kalman. 
  #-------------------------------------------------------------#
  if(!is.null(varLabels$conv)){
    
    ts_list <- setupConvolve(
      ts_list       = ts_list, 
      varLabels     = varLabels, 
      conv_length   = ctrlOpts$conv_length, 
      conv_interval = ctrlOpts$conv_interval
    )
    
  }
  #-------------------------------------------------------------#
  
  
  #-------------------------------------------------------------#
  # (2) standardize (including conv_vars post convolution)
  #-------------------------------------------------------------#
  # If neccessary, standardize variables at t.
  #   * conv_vars at t should be standardized
  #   * categorical variables should not be standardized
  #       Note: conv_vars should NOT be labeled as categorical
  #             add check for this and throw error. Also, add
  #             note to exogenous documentation that categorical
  #             vars should be specified in categorical arg.
  #   * bilinear vars should not be standardized.
  #-------------------------------------------------------------# 
  if(ctrlOpts$standardize){
    
    vars_to_scale <- setdiff(varLabels$coln, c(varLabels$mult, varLabels$lagg, varLabels$categorical))
    
    ts_list <- lapply(ts_list, function(df) {
      
      df[,vars_to_scale] <- scale(df[,vars_to_scale])
      df
      
    })
  
  }
  #-------------------------------------------------------------# 
  
  
  #-------------------------------------------------------------#
  # (3) lag data
  #-------------------------------------------------------------#
  #  * by default, lagged conv_vars are created, 
  #      and they must be subsequently removed.
  #-------------------------------------------------------------#
  ts_list <- lapply(ts_list, function(df){
    
    first           <- df[1:(nrow(df)-1), ] 
    second          <- df[2:(nrow(df)  ), ]
    ts_lc           <- data.frame(first, second)
    colnames(ts_lc) <- c(paste0(colnames(df), "lag"), colnames(df))
    
    # if neccessary, remove lagged conv_vars
    all_poss_lag_conv_vars <- paste0(setdiff(varLabels$conv, varLabels$lagg),"lag")
    non_spec_lag_conv_vars <- setdiff(all_poss_lag_conv_vars, varLabels$conv)
    
    
    if(length(non_spec_lag_conv_vars) > 0){
      
      ts_lc <- ts_lc[,!(colnames(ts_lc) %in% non_spec_lag_conv_vars), drop = FALSE]
      
    }
    
    ts_lc
    
  })
  #-------------------------------------------------------------#
  
  
  
  #-------------------------------------------------------------#
  # (4) create bilinear variables
  #-------------------------------------------------------------#
  #  
  #-------------------------------------------------------------#
  if(!is.null(varLabels$mult)){
    
    mult_vars_list <- strsplit(varLabels$mult, "_by_", fixed = TRUE)
    
    new_mult_df <- function(v1,v2, data){
      
      df <- data.frame(data[,v1]*data[,v2])
      colnames(df) <- paste0(v1,"_by_", v2)
      
      if(ctrlOpts$mean_center_mult) {
        df[] <- scale(df[,1], center = TRUE, scale = FALSE)
      } 
      
      df
    }
    
    ts_list <- lapply(ts_list, function(df){
      
      cbind(df, do.call("cbind", 
        lapply(mult_vars_list, function(m){new_mult_df(m[1],m[2], df)})
      ))
      
    })
  }
  #-------------------------------------------------------------# 

  

  #-------------------------------------------------------------#
  # Create directories and catch errors. 
  #-------------------------------------------------------------#
  if (is.null(ctrlOpts$out)){
    
    cat(
      "gimme MESSAGE: No output directory specified. ",
      "All output should be directed to an object.", "\n"
    )
    
  } else if (file_test(op = "-d", ctrlOpts$out) == FALSE) {
    
    cat(
      "gimme MESSAGE: specified output directory doesn't exist. ",
      "Attempting to create now.", "\n"
      
    )
    dir.create(ctrlOpts$out, recursive = TRUE)
    
    if (dir.exists(ctrlOpts$out)){
      
      cat(
        "gimme MESSAGE: output directory successfully created.", "\n"
      )
      
    } else {
      
      stop(
        "gimme ERROR: unable to create output directory. ",
        "Please specify a different file path"
      )
      
    }
    
  }
  
  individual        <- file.path(ctrlOpts$out, "individual")
  
  subgroup_dir      <- file.path(ctrlOpts$out, "subgroup")
  
  if (ctrlOpts$subgroup & !is.null(ctrlOpts$out))  {
    
    dir.create(subgroup_dir, showWarnings = FALSE)
    
  }
  
  if (!ctrlOpts$agg & !is.null(ctrlOpts$out)) {
    
    dir.create(individual, showWarnings = FALSE)
    
  }
  
  #-------------------------------------------------------------#
  
  
  
  #-------------------------------------------------------------#
  # Final data checks
  #-------------------------------------------------------------#
  
  # todo
  # reorder colnames? 
  # throw warning if varLabels$coln does not equal colnames(ts_list[[1]])
  ts_list <- lapply(ts_list, function(x){x[,varLabels$coln]})
   
  
  #-------------------------------------------------------------#
  # More data checks.
  #-------------------------------------------------------------#


  n_subjects   <- length(ts_list)
  cols         <- numeric()
  missingCols  <- numeric()
  constantCols <- logical()
  numericCols  <- logical()
  
  # check for obvious errors in data
  for (k in 1:length(ts_list)){
    data.file <- ts_list[[k]]
    cols[k]   <- ncol(data.file)
    missingCols[k] <- sum(colSums(is.na(data.file)) < nrow(data.file))
    constantCols[k] <- any(apply(data.file, 2, sd, na.rm = TRUE) == 0)
    numericCols[k]  <- any(apply(data.file, 2, is.numeric) == FALSE)
  }
  
  
  if (n_subjects != 1) {
    if (sd(cols) != 0) {
      stop(paste0('gimme ERROR: not all data files have the same number of columns. ',
                  'Please fix or remove file before continuing.'))
    }
    if (sd(missingCols) != 0) {
      stop(paste0('gimme ERROR: at least one data file contains a column with all NA. ',
                  'Please fix or remove files listed below before continuing. \n', 
                  paste0(names(ts_list)[missingCols != cols], collapse = "\n")))
    }
    if (any(cols != missingCols)) {
      stop(paste0('gimme ERROR: at least one data file contains a column with all NA. ',
                  'Please fix or remove file before continuing.'))
    }  
    if (any(constantCols == TRUE)){
      stop(paste0('gimme ERROR: at least one data file contains a column with constant values. ',
                  'Please fix or remove files listed below before continuing. \n', 
                  paste0(names(ts_list)[constantCols == TRUE], collapse = "\n")))
    }
    if (any(numericCols == TRUE)){
      stop(paste0('gimme ERROR: at least one data file contains a column with non-numeric values. ',
                  'Please fix or remove files listed below before continuing. \n', 
                  paste0(names(ts_list)[numericCols == TRUE], collapse = "\n")))
    }
  } 
  if (n_subjects == 1 & !ctrlOpts$ind) {
    stop(paste0('gimme ERROR: only one subject detected in data directory. ',
                'Please use indSEM function instead.'))
  }
  
  
  return(ts_list)
  
}