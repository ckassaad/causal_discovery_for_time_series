#' Do some preliminary checks on the data.
#' @keywords internal
setupPrelimDataChecks <- function(df){
  
  data <- lapply(df, function(dfi){
    first           <- dfi[1:(nrow(dfi)-1), ] 
    second          <- dfi[2:(nrow(dfi)  ), ]
    data            <- data.frame(first, second)
    colnames(data)  <- c(paste0(colnames(dfi), "lag"), colnames(dfi))
    data
  })


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

  
  return(df)
  
}