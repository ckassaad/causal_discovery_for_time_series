#' @keywords internal
setup <- function (data,
                   sep,
                   header,
                   out,
                   plot,
                   ar,
                   paths,
                   exogenous,
                   mult_vars,
                   mean_center_mult,
                   standardize,
                   subgroup,
                   agg,
                   ind,
                   groupcutoff,
                   subcutoff, 
                   conv_vars, 
                   conv_length, 
                   conv_interval,
                   lv_model = NULL,
                   lv_scores = "regression",       
                   lv_estimator = "miiv",            
                   lv_miiv_scaling  = "first.indicator", 
                   lv_final_estimator = "miiv") {
  

    
  #-------------------------------------------------------------#
  # NOTES and TODOS
  #-------------------------------------------------------------#
  #
  # Add error if any variables contain the following strings:
  #   "_by_".
  # Add error if ex-lag is specified, tell user function is 
  #   deprecated. 
  # Make print method for ctrlOpts class, similar to varLabels.
  #-------------------------------------------------------------#
  
  
  
  #-------------------------------------------------------------#
  # Assemble static control options.
  #
  #  There is an associated print method for the ctrlOpts
  #    class. For example, print(ctrlOpts). [todo]
  #
  #-------------------------------------------------------------#

  
  ctrlOpts <- list(
    sep              = sep,
    header           = header,
    out              = out,
    ind_dir          = file.path(out, "individual"),
    sub_dir          = file.path(out, "subgroup"),
    plot             = plot,
    ar               = ar,
    mean_center_mult = mean_center_mult,
    standardize      = standardize,
    subgroup         = subgroup,
    agg              = agg,
    ind              = ind,
    groupcutoff      = groupcutoff,
    subcutoff        = subcutoff, 
    conv_length      = conv_length, 
    conv_interval    = conv_interval,
    diagnostics      = NULL,
    paths            = paths
  )
  
  class(ctrlOpts) <- "ctrlOpts"
  
  # if (diagnostics) { print(ctrlOpts) }
  
  #-------------------------------------------------------------#
  
  
  #-------------------------------------------------------------#
  # Create list of dataframes.
  #   # Notes:
  #   #  Now all datasets share the same column order.
  #   #
  #   # Do this here:
  #   # remove lagged conv_vars (if neccessary)
  #   # remove lagged exog_vars (if neccessary)
  #
  #-------------------------------------------------------------#
  
  ts_list <- setupDataLists(data = data, ctrlOpts = ctrlOpts, lv_model = lv_model)
  
  
  #-------------------------------------------------------------#
  # lv_gimme
  #-------------------------------------------------------------#

  if(is.null(lv_model)){
    
     lvgimme <- NULL
    
  } else {
    
    ts_list <- setupPrelimDataChecks(ts_list)
    
    ts_list_obs <- ts_list
    
    if(length(lv_model) == 1){
      
       writeLines(paste0("lvgimme: measurement model invariant across subjects"))
       lv_model_all <- replicate(length(ts_list_obs), lv_model)
       same_meas_model  <- TRUE
       
    } else if (length(lv_model) == length(ts_list_obs)){
      
      writeLines(paste0("lvgimme: measurement model varies across subjects"))
      
      if(lv_miiv_scaling == "group"){
        stop(paste0("error: group scaling not possible with individual-level measurement models"))
      }
      
      lv_model_all <- lv_model
      
    } else {
      
      stop(paste0("error: lv_model must be of length 1 or ",length(ts_list_obs)))
      
    }
    
    # zf: uncomment once changes pushed through final MIIV estimation
    # first let's look for single indicator LVs and extract the raw data:
    single_indicator_lvs <- lapply(seq_along(lv_model_all), function(i){
      
      pt  <- lavaan::lavParTable(lv_model_all[[i]])
      lvs <- unique(pt[pt$op=="=~","lhs"])
      
      lv_list <- list(); cnt <- 1
               
       for(l in 1:length(lvs)){
          if(length(pt[pt$op=="=~" & pt$lhs == lvs[l],"rhs"]) == 1){
            
            keep <- pt[pt$op=="=~" & pt$lhs == lvs[l],"rhs"]
            tmp <- ts_list_obs[[i]][,keep, drop = F]
            colnames(tmp) <- lvs[l]
            lv_list[[cnt]] <- tmp
            cnt <- cnt + 1
          } 
       }        
      
       do.call( "cbind", lv_list)

    })
    
    
    lv_names <- list()
    ov_names <- list()
    
    for(j in 1:length(lv_model_all)){
      
      pt  <- lavaan::lavParTable(lv_model_all[[j]])
      lvs <- unique(pt[pt$op=="=~","lhs"])
      
      for(k in 1:length(lvs)){
        
        if(length(pt[pt$op=="=~" & pt$lhs == lvs[k],"rhs"]) == 1){
          
          lv_names_v <- lvs[k]
          ov_names_v <- pt[pt$op=="=~" & pt$lhs == lvs[k],"rhs"]
          
          if(length(ov_names_v) > 1){
            stop(paste0("gimme: internal error processing variable names."))
          }
        
          colnames(ts_list_obs[[j]]) <- gsub(
            ov_names_v, 
            lv_names_v, 
            colnames(ts_list_obs[[j]])
          )
          
        }
      }
    }

    
    
    # # since these are observed variables every model should have
    # # the same set. otherwise this will cause an error later. 
    
    # now we need to remove the single indicator LVs from the model
    lv_model_ind <- lapply(lv_model_all, function(x){
      
      pt       <- lavaan::lavParTable(x)
      lvs      <- unique(pt[pt$op=="=~","lhs"])
      mod.list <- list()
      cnt      <- 1L
      
      for(j in 1:length(lvs)){
        
        if (length(pt[pt$op=="=~" & pt$lhs == lvs[j],"rhs"]) == 1){
          
          # remove this error once the functionality is pushed through the final MIIV estimation.
          #stop(paste0("gimme ERROR: factors with only one indicator not currently supported."))
          
        } else if (length(pt[pt$op=="=~" & pt$lhs == lvs[j],"rhs"]) == 2){
          
           stop(paste0("gimme ERROR: factors with only two indicators not currently supported."))
          
        } else {
        
          mod.list[[cnt]] <- paste0(paste0(lvs[j],"=~"), paste0(pt[pt$op=="=~" & pt$lhs == lvs[j],"rhs"],collapse="+"),collapse="")
          cnt <- cnt + 1
          
        }
        
      }
      
      mod.list
      
    })
    

    model_list_ptech <- getModelList(
      data = ts_list_obs, 
      model = lv_model_ind, 
      scaling = lv_miiv_scaling, 
      lag = FALSE
    )
    
    model_list_dfa <- getModelList(
      data = ts_list_obs, 
      model = lv_model_ind, 
      scaling = lv_miiv_scaling, 
      lag = TRUE
    )

    fs_info <- lapply(seq_along(ts_list_obs), function(i){
      
      lapply(seq_along(model_list_ptech[[i]]), function(j){
        factorScores(
          data  = ts_list_obs[[i]], 
          model = model_list_ptech[[i]][j], 
          lv_scores = lv_scores, 
          lv_estimator = lv_estimator 
        )   
      })
      
      
    })
  
    ts_list <- lapply(seq_along(fs_info), function(i){
      as.data.frame(do.call("cbind", lapply(fs_info[[i]],"[[", "fs")))
    })   
    
    # add back in any single indicator lvs
    # zf: uncomment once changes pushed through final MIIV estimation
    if(!all(unlist(lapply(single_indicator_lvs,function(x){is.null(x)})))){
      ts_list <- lapply(seq_along(ts_list), function(i){
        cbind(ts_list[[i]], single_indicator_lvs[[i]])
      })
    }

    names(ts_list) <- names(ts_list_obs)
    
    lvgimme <- list(
      ts_list_obs = ts_list_obs,
      lv_model_ind = lv_model_ind,
      model_list_ptech = model_list_ptech,
      model_list_dfa = model_list_dfa,
      fs_info = fs_info
    )
         
  }
    

  #-------------------------------------------------------------#
  # Assemble different variable sets.
  #
  #  There is an associated print method for the varLabels
  #    class. For example, print(varLabels).
  #
  #-------------------------------------------------------------#
  
  ### Distinguish between lagged and contemporaenous exogenous variables
  exog_con <- exogenous[regexpr("&lag", exogenous)<0]
  exog_lag <- sub("&lag", "", exogenous[regexpr("&lag", exogenous)>0])
  exogenous <- c(exog_lag,exog_con)
  
  orig <- colnames(ts_list[[1]])
  uexo <- unique(exogenous)
  conv <- conv_vars
  lagg <- paste0(setdiff(orig,unique(exog_con, conv)), "lag")
  mult <- setupMultVarNames(mult_vars)
  exog <- unique(c(uexo, mult, lagg))
  endo <- setdiff(orig, exog) # only true if ar = TRUE
  catg <- NULL
  stnd <- if(standardize) setdiff(c(endo,exog), c(catg, conv_vars)) else NULL
  #coln <- c(endo,exog) # future column names of data
  coln <- unique(c(lagg, orig, mult))
  

  varLabels <- list(
    orig = orig,
    uexo = uexo,
    lagg = lagg,
    conv = conv,
    mult = mult,
    exog = exog,
    endo = endo,
    catg = catg,
    stnd = stnd,
    coln = coln,
    exog_lag = exog_lag,
    exog_con = exog_con
  )
  
  class(varLabels) <- "varLabels"
  
  #if (diagnostics) { print(varLabels) }
  
 
  
  #-------------------------------------------------------------#
  
  
  #-------------------------------------------------------------#
  # Final data manipulation.
  #-------------------------------------------------------------#
  # (1) early data checks
  # (2) convolve contemporaneous 
  # (3) standardize (including conv_vars post convolution)
  # (4) lag data (creating any lagged conv_vars)
  # (5) remove lagged conv_vars (if neccessary)
  # (6) create bilinear (may require lagged variables)
  # (7) late data checks
  #-------------------------------------------------------------#
  
  ts_list <- setupTransformData(
    ts_list   = ts_list,
    varLabels = varLabels,
    ctrlOpts  = ctrlOpts
  )

  #-------------------------------------------------------------#
  
  #-------------------------------------------------------------#
  # Plotting information. Why is this here?
  #-------------------------------------------------------------#
  #
  # this needs to be fixed as n_lagged is meaningless.
  # if (plot) {
  #   plot.names <- varnames[(n_lagged+1):(n_vars_total)] 
  # } else {
  #   plot.names <- ""
  # }
  #
  #-------------------------------------------------------------#
  
  
  #-------------------------------------------------------------#
  # Prepare paths argument if semigimme is specified
  #-------------------------------------------------------------#
  if (!is.null(ctrlOpts$paths)) {

    pathInfo <- setupPrepPaths(ctrlOpts$paths, varLabels, ctrlOpts)
    remove   <- pathInfo$remove  
    paths    <- pathInfo$paths
    
  } else {
    
    paths    <- ctrlOpts$paths
    remove   <- NULL
    
  }
  #-------------------------------------------------------------#
  
  
  #-------------------------------------------------------------#
  # Prepare paths argument if semigimme is specified
  #-------------------------------------------------------------#
  pathList <- setupBaseSyntax(paths, varLabels, ctrlOpts)
  
  #-------------------------------------------------------------#
  
  
  dat <- list("ar" = ctrlOpts$ar, 
              "out"= ctrlOpts$out,
              "plot" = ctrlOpts$plot,
              "subgroup" = ctrlOpts$subgroup,
              "agg" = ctrlOpts$agg,
              "n_subj" = length(ts_list),
              "n_lagged" = length(varLabels$lagg),
              "n_exog_lag" = length(exog_lag),
              "n_exog" = length(varLabels$uexo),
              "n_bilinear" = length(varLabels$mult),
              "n_endog"  = length(varLabels$endo),
              "n_exog_total" = length(c(varLabels$uexo, varLabels$mult)),
              "n_vars_total" = length(varLabels$coln),
              "n_contemporaneous" = length(varLabels$coln)-length(varLabels$lagg),
              "nonsense_paths" = pathList$nonsense.paths,
              "varnames" = varLabels$coln,
              "lvarnames" = varLabels$coln,
              "cutoffind" = qchisq(.99, 1),
              "subgroup_dir" = ctrlOpts$sub_dir,
              "ind_dir"   = ctrlOpts$ind_dir,
              "syntax"     = pathList$syntax,
              "candidate_paths" = pathList$candidate.paths,
              "candidate_corr" = pathList$candidate.corr,
              "fixed_paths"  = pathList$fixed.paths,
              "group_cutoff" = ctrlOpts$groupcutoff,
              "sub_cutoff" = ctrlOpts$subcutoff,
              "ts_list" = ts_list,
              "standardize" = ctrlOpts$standardize,
              "file_order" =  data.frame(
                index = c(seq(1:length(ts_list))), 
                names = c(names(ts_list)), 
                stringsAsFactors = FALSE
              ),
              "chisq_cutoff_mi_epc" = qchisq(1-(.05/((length(varLabels$coln)*(length(varLabels$coln)-1)/2)*length(ts_list))), 1),
              "varLabels" = varLabels,
              "ctrlOpts"  = ctrlOpts,
              "lvgimme"   = lvgimme
    )
  return(dat)
}