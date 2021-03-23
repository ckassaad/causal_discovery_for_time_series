#' Create a list of dataframes
#'
#' @param data a list or directory.
#' @param ctrlOpts a lit of control options.
#' @keywords internal
setupDataLists <- function(data, ctrlOpts  = NULL, lv_model = NULL){
  
  ## throw error if data argument is null
  if (is.null(data)){
    stop(paste0(
      "gimme ERROR: neither a data directory nor a data list is specified. ",
      "Please either specify a directory of data or a list of individual data files."
    )
    )
  }
  
  
  #-------------------------------------------------------------#
  # If the data is already in list form.
  #-------------------------------------------------------------#
  
  if (is.list(data)){
  
    ts_list  <- list()
  
    # if the user-supplied list does not have names, add names 
    if(is.null(names(data))){ 
      
      names(data) <- paste0("subj", 1:length(data)) 
      
    }
  
    ts_list  <- data

  #-------------------------------------------------------------#
  # If the data is not a list (we assume it is a directory).
  #-------------------------------------------------------------#  
  } else if (!is.list(data)){
  
    files <- list.files(data, full.names = TRUE) 
    
    #-------------------------------------------------------------#
    # Throw errors specific to data specified as a directory.
    #-------------------------------------------------------------# 
    if (is.null(ctrlOpts$sep)){
      stop(paste0(
        "gimme ERROR: a data directory is specified but a ctrlOpts$sep argument is not. ",
        "Please specify a ctrlOpts$sep argument before continuing."
      ))
    }
    
    
    if (is.null(ctrlOpts$header)){
      stop(paste0(
        "gimme ERROR: a data directory is specified but a ctrlOpts$header argument is not. ",
        "Please specify a logical value for ctrlOpts$header before continuing."
      ))
    }
  
    #-------------------------------------------------------------#
    # Create a list of dataframes.
    #-------------------------------------------------------------# 
    
    ts_list <- list()
    for (i in 1:length(files)){
      ts_list[[i]] <- read.table(files[i], sep = ctrlOpts$sep, header = ctrlOpts$header)
    }
    
    names(ts_list) <- tools::file_path_sans_ext(basename(files))
  
  } else {
  
    stop(paste0(
      "gimme ERROR: Format of data argument not recognized. "
    ))

  }
  
  #-------------------------------------------------------------#
  # Ensure all datafiles share the same column order.
  #
  #   # Now we also make all datasets have the same column order.
  #-------------------------------------------------------------# 
  varnames       <- colnames(ts_list[[1]])
  n_orig_vars    <- ncol(ts_list[[1]])
  
  if(is.null(varnames)){
    varnames <- c(paste0("V", seq(1,n_orig_vars)))
    ts_list <- lapply(ts_list, function(x) { colnames(x)<-varnames;  x })
  } else {
    if(is.null(lv_model)){
      ts_list <- lapply(ts_list, function(x) { x[,varnames]})
    }
  }
  
    
  return(ts_list)
  
}