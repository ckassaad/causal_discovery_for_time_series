#' Do some final checks on the data.
#' @keywords internal
setupFinalDataChecks <- function(data){
  # satisfy CRAN checks
  ts_list = NULL
  ind = NULL
  
  n_subjects   <- length(data)
  cols         <- numeric()
  missingCols  <- numeric()
  constantCols <- logical()
  numericCols  <- logical()
  
  # check for obvious errors in data
  for (k in 1:length(data)){
    data.file <- data[[k]]
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
  if (n_subjects == 1 & !ind) {
    stop(paste0('gimme ERROR: only one subject detected in data directory. ',
                'Please use indSEM function instead.'))
  }
  
  return(data)
  
}