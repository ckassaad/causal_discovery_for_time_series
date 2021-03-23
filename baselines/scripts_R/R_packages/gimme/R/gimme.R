#' @name gimmeSEM
#' @aliases gimme gimmeSEM
#' @title Group iterative multiple model estimation.
#' @description This function identifies structural equation models for each
#' individual that consist of both group-level and individual-level paths.
#' @usage
#' gimmeSEM(data        = NULL,
#'          out         = NULL,
#'          sep         = NULL,
#'          header      = NULL,
#'          ar          = TRUE,
#'          plot        = TRUE,
#'          subgroup    = FALSE,
#'          sub_feature = "lag & contemp",
#'          sub_method = "Walktrap",
#'          sub_sim_thresh    = "lowest", 
#'          confirm_subgroup = NULL,
#'          paths       = NULL,
#'          exogenous = NULL,
#'          conv_vars   = NULL,
#'          conv_length = 16, 
#'          conv_interval = 1,
#'          mult_vars   = NULL,
#'          mean_center_mult = FALSE,
#'          standardize = FALSE,
#'          groupcutoff = .75,
#'          subcutoff   = .5,
#'          diagnos     = FALSE, 
#'          ms_allow         = FALSE,
#'          ms_tol           = 1e-5,
#'          lv_model         = NULL, 
#'          lv_estimator     = "miiv",     
#'          lv_scores        = "regression",       
#'          lv_miiv_scaling  = "first.indicator", 
#'          lv_final_estimator = "miiv",
#'          lasso_model_crit    = NULL, 
#'          hybrid = FALSE,
#'          dir_prop_cutoff =0)
#' @param data The path to the directory where the data files are located,
#' or the name of the list containing each individual's time series. Each file
#' or matrix must contain one matrix for each individual containing a T (time)
#' by p (number of variables) matrix where the columns represent variables and
#' the rows represent time. Individuals must have the same variables (p)
#' but can have different lengths of observations (T).
#' @param out The path to the directory where the results will be stored
#' (optional). If specified,
#' a copy of output files will be replaced in directory. If directory at
#' specified path does not exist, it will be created.
#' @param sep The spacing of the data files. 
#' "" indicates space-delimited,
#' "/t" indicates tab-delimited, "," indicates comma delimited. Only necessary
#' to specify if reading data in from physical directory.
#' @param header Logical. Indicate TRUE for data files with a header. Only
#' necessary to specify if reading data in from physical directory.
#' @param ar Logical. If TRUE, begins search for group model with
#' autoregressive (AR) paths freed for estimation. If ms_allow=TRUE, it is recommended
#' to set ar=FALSE.  Multiple solutions are unlikely to be found when ar=TRUE.  Defaults to TRUE.
#' @param paths \code{lavaan}-style syntax containing paths with which
#' to begin model estimation (optional). That is, Y~X indicates that Y
#' is regressed on X, or X predicts Y. If no header is used,
#' then variables should be referred to with V followed (with no separation)
#' by the column number. If a
#' header is used, variables should be referred to using variable names.
#' To reference lag variables, "lag" should be added to the end of the variable
#' name with no separation. Defaults to NULL. 
#' @param exogenous Vector of variable names to be treated as exogenous (optional).
#' That is, exogenous variable X can predict Y but cannot be predicted by Y.
#' If no header is used, then variables should be referred to with V followed
#' (with no separation) by the column number.  If a header is used, variables should be referred 
#' to using variable names. The default for exogenous variables is that lagged effects of the exogenous 
#' variables are not included in the model search.  If lagged paths are wanted, "&lag" should be added to the end of the variable
#' name with no separation. Defaults to NULL.
#' @param conv_vars Vector of variable names to be convolved via smoothed Finite Impulse 
#' Response (sFIR). Note, conv_vars are not not automatically considered exogenous variables.
#' To treat conv_vars as exogenous use the exogenous argument. Variables listed in conv_vars 
#' must be binary variables. If there is missing data in the endogenous variables their values 
#' will be imputed for the convolution operation only. Defaults to NULL. ### If there are multiple 
#' variables listed in conv_vars they are not used in the convolution of additional conv_vars.## 
#' You can't do lagged variables.
#' @param conv_length Expected response length in seconds. For functional MRI BOLD, 16 seconds (default) is typical
#' for the hemodynamic response function. 
#' @param conv_interval Interval between data acquisition. Currently must be a constant. For 
#' fMRI studies, this is the repetition time. Defaults to 1. 
#' @param mult_vars Vector of variable names to be multiplied to explore bilinear/modulatory
#' effects (optional). All multiplied variables will be treated as exogenous (X can predict
#' Y but cannot be predicted by Y). Within the vector, multiplication of two variables should be
#' indicated with an asterik (e.g. V1*V2). If no header is used, variables should be referred to with 
#' V followed by the column number (with no separation). If a header is used, each variable should be
#' referred to using variable names. If multiplication with the lag 1 of a variable is desired, the 
#' variable name should be followed by "lag" with no separation (e.g. V1*V2lag). 
#' @param mean_center_mult Logical. If TRUE, the variables indicated in mult_vars will be mean-centered
#' before being multiplied together. Defaults to FALSE. 
#' @param standardize Logical. If TRUE, all variables will be standardized to have a mean of zero and a
#' standard deviation of one. Defaults to FALSE. 
#' @param plot Logical. If TRUE, graphs depicting relations among variables
#' of interest will automatically be
#' created. Solid lines represent contemporaneous relations (lag 0) and dashed lines reflect 
#' lagged relations (lag 1). For individual-level plots, red paths represent positive weights
#' and blue paths represent negative weights. Width of paths corresponds to estimated path weight.
#' For the group-level plot, black represents group-level paths, grey represents
#' individual-level paths, and (if subgroup = TRUE)
#' green represents subgroup-level paths. For the group-level plot,
#' the width of the edge corresponds to the count. Defaults to TRUE.
#' @param subgroup Logical. If TRUE, subgroups are generated based on
#' similarities in model features using the \code{walktrap.community}
#' function from the \code{igraph} package. When ms_allow=TRUE, subgroup
#' should be set to FALSE.  Defaults to FALSE. 
#' @param confirm_subgroup Dataframe. Option only available when subgroup = TRUE. Dataframe should contain two columns. The first
#' column should specify file labels (the name of the data files without file extension), 
#' and the second should contain integer values (beginning at 1) 
#' specifying the subgroup membership for each individual.
#' function from the \code{igraph} package. Defaults to TRUE. 
#' @param sub_feature Option to indicate feature(s) used to subgroup individuals. Defaults to
#' "lag & contemp" for lagged and contemporaneous, which is the original method. Can use 
#' "lagged" or "contemp" to subgroup solely on features related to lagged and contemporaneous 
#' relations, respectively.
#' @param sub_method Community detection method used to cluster individuals into subgroups. Options align 
#' with those available in the igraph package: "Walktrap" (default), "Infomap", "Louvain", "Edge Betweenness", 
#' "Label Prop", "Fast Greedy", "Leading Eigen", and "Spinglass". 
#' @param sub_sim_thresh Threshold for inducing sparsity in similarity matrix. Options are: the percent of edges 
#' in the similarity matrix to set to zero (e.g., .25 would set the lower quartile), "lowest" (default) subtracts 
#' the minimum value from all values, and "search" searches across thresholds to arrive at one providing highest modularity.  
#' @param groupcutoff Cutoff value for group-level paths. Defaults to .75,
#' indicating that a path must be significant across 75\% of individuals to be
#' included as a group-level path.
#' @param subcutoff Cutoff value for subgroup- level paths. Defaults to .5,
#' indicating that a path must be significant across at least 50\% of the
#' individuals in a subgroup to be considered a subgroup-level path.
#' @param diagnos Logical.If TRUE provides internal output for diagnostic purposes. Defaults to FALSE. 
#' @param ms_allow Logical. If TRUE provides multiple solutions when more than one path has identical 
#' modification index values.  When ms_allow=TRUE, it is recommended
#' to set ar=FALSE.  Multiple solutions are unlikely to be found when ar=TRUE.  Additionally,
#' subgroup should be set to FALSE.  Output files for individuals with multiple solutions will represent the last solution 
#' found for the individual, not necessarily the best solution for the individual.
#' @param ms_tol Precision used when evaluating similarity of modification indices when ms_allow = TRUE.  We recommend
#' that ms_tol not be greater than the default, especially when standardize=TRUE.     
#' Defaults to 1e-5.
#' @param lv_model Invoke latent variable modeling by providing the measurement model syntax here. lavaan
#' conventions are used for relating observed variables to factors. Defaults to NULL.
#' @param lv_estimator Estimator used for factor analysis. Options are "miiv" (default), "pml" (pseudo-ML) or "svd".
#' @param lv_scores Method used for estimating latent variable scores from parameters obtained from the factor analysis 
#' when lv_model is not NULL. Options are: "regression" (Default), "bartlett".
#' @param lv_miiv_scaling Type of scaling indicator to use when "miiv" selected for lv_estimator. Options are
#' "first.indicator" (Default; the first observed variable in the measurement equation is used), "group" 
#' (best one for the group), or "individual" (each individual has the best one for them according to R2). 
#' @param lv_final_estimator Estimator for final estimations. "miiv" (Default) or "pml" (pseudo-ML). 
#' @param lasso_model_crit When not null, invokes multiLASSO approach for the GIMME model search procedure. Arguments 
#' indicate the model selection criterion to use for model selection: 'bic' (select on BIC), 'aic', 'aicc', 'hqc', 'cv' (cross-validation). 
#' @param hybrid Logical. If TRUE, enables hybrid-VAR models where both directed contemporaneous paths and contemporaneous 	
#' covariances among residuals are candidate relations in the search space. Defaults to FALSE.
#' @param dir_prop_cutoff Option to require that the directionality of a relation has to be higher than the reverse direction for a prespecified proportion of indivdiuals.  
#' @details
#'  In main output directory:
#'  \itemize{
#'  \item{\strong{indivPathEstimates}} {Contains estimate, standard error,
#'  p-value, and z-value for each path for each individual.
#'  If subgroup = TRUE and subgroups are found, then a column is present
#'  containing the subgroup membership for each individual. Also contains the
#'  level at which each path was estimated: group, subgroup, or individual.}
#'  \item{\strong{summaryFit}} {Contains model fit information for individual-
#'  level models. If subgroups are requested, this file also contains the
#'  subgroup membership for each individual.}
#'  \item{\strong{summaryPathCountMatrix}} Contains counts of total number of
#'  paths, both contemporaneous and lagged, estimated for the sample. The row
#'  variable is the outcome and the column variable is the predictor variable.
#'  \item{\strong{summaryPathCounts}} {Contains summary count information for
#'  paths identified at the group-, subgroup (if subgroup = TRUE), and
#'  individual-level.}
#'  \item{\strong{summaryPathsPlot}} {Produced if plot = TRUE. Contains figure
#'  with group, subgroup (if subgroup = TRUE), and individual-level paths
#'  for the sample. Black paths are group-level, green paths are subgroup-level,
#'  and grey paths are individual-level, where the thickness of the line
#'  represents the count.}
#'  }
#'  In subgroup output directory (if subgroup = TRUE):
#'  \itemize{
#'  \item{\strong{subgroup\emph{k}PathCounts}} Contains counts of relations
#'  among lagged and contemporaneous variables for the \strong{\emph{k}}th
#'  subgroup.
#'  \item{\strong{subgroup\emph{k}Plot}} Contains plot of group, subgroup,
#'  and individual level paths for the \strong{\emph{k}}th subgroup.
#'  Black represents group-level paths, grey represents individual-level paths,
#'  and green represents subgroup-level paths.
#'  }
#'  Note: if a subgroup of size n = 1 is discovered, subgroup-level output is
#'  not produced. \cr
#'  In individual output directory (where \strong{\emph{id}} represents the
#'  original file name for each individual):
#'  \itemize{
#'  \item{\strong{\emph{id}Betas}} Contains individual-level estimates
#'   of each path for each individual.
#'  \item{\strong{\emph{id}StdErrors}} Contains individual-level standard errors
#'  for each path for each individual.
#'  \item{\strong{\emph{id}Plot}} Contains individual-level plots. Red paths
#'  represent positive weights and blue paths represent negative weights.
#' }
#' @references Gates, K.M. & Molenaar, P.C.M. (2012). Group search algorithm
#' recovers effective connectivity maps for individuals
#' in homogeneous and heterogeneous samples. NeuroImage, 63, 310-319.
#' @references Lane, S.T. & Gates, K.M. (2017). Automated selection of robust
#' individual-level structural equation models for time series data.
#' Structural Equation Modeling.
#' @references Adriene M. Beltz & Peter C. M. Molenaar (2016) Dealing 
#' with Multiple Solutions in Structural Vector Autoregressive Models, 
#' Multivariate Behavioral Research, 51:2-3, 357-373.
#' @author Stephanie Lane and Zachary Fisher
#' @examples
#'  \dontrun{
#' paths <- 'V2 ~ V1
#'           V3 ~ V4lag'
#'
#' fit <- gimmeSEM(data     = simData,
#'                 out      = "C:/simData_out",
#'                 subgroup = TRUE, 
#'                 paths    = paths)
#'
#' print(fit, mean = TRUE)
#' print(fit, subgroup = 1, mean = TRUE)
#' print(fit, file = "group_1_1", estimates = TRUE)
#' print(fit, subgroup = 2, fitMeasures = TRUE)
#' plot(fit, file = "group_1_1")
#'  }
#' @keywords gimmeSEM
#' @export gimme gimmeSEM

gimmeSEM <- gimme <- function(data             = NULL,
                              out              = NULL,
                              sep              = NULL,
                              header           = NULL,
                              ar               = TRUE,
                              plot             = TRUE,
                              subgroup         = FALSE,
                              sub_feature      = "lag & contemp",
                              sub_method       = "Walktrap",
                              sub_sim_thresh   = "lowest", 
                              confirm_subgroup = NULL,
                              paths            = NULL,
                              exogenous        = NULL,
                              conv_vars        = NULL,
                              conv_length      = 16, 
                              conv_interval    = 1, 
                              mult_vars        = NULL,
                              mean_center_mult = FALSE,
                              standardize      = FALSE,
                              groupcutoff      = .75,
                              subcutoff        = .5,
                              diagnos          = FALSE,
                              ms_allow         = FALSE,
                              ms_tol           = 1e-5,
                              lv_model         = NULL, 
                              lv_estimator     = "miiv",             # c("miiv", "pml")
                              lv_scores        = "regression",       # c("regression", "bartlett")
                              lv_miiv_scaling  = "first.indicator",  # c("group", "individual")
                              lv_final_estimator = "miiv",
                              lasso_model_crit = NULL, 
                              hybrid = FALSE,
                              dir_prop_cutoff = 0){          # c("miiv", "pml")
  
  # satisfy CRAN checks
  ind     = NULL
  varnames = NULL
  lvarnames = NULL
  sub_membership = NULL

  setupConvolve = NULL
  ts = NULL
  setupFinalDataChecks = NULL
  
  
  #Error check for ms_allow
  if(ms_allow & subgroup){
      stop(paste0("gimme ERROR: Subgrouping is not available for ms-gimme.",
                  " Please ensure that subgroup=FALSE if ms_allow=TRUE"))
  }
  
  if(ms_allow & ar){
    writeLines("gimme WARNING: Multiple solutions are not likely when ar=TRUE.",
                " We recommend setting ar to FALSE if using ms_allow.")
  }
  
  #Error check for hybrid
  if(hybrid & !ar){
    stop(paste0("gimme ERROR: Autoregressive paths have to be open for hybrid-gimme.",
                " Please ensure that ar=TRUE if hybrid=TRUE."))
  }
  
   sub_membership = NULL
   
   # if !is.null(lasso_model_crit)
   # {
   #   if(!is.null(mult_vars)){
   #   ml <- strsplit(mult_vars, "*", fixed = TRUE)
   #   
   #   # identify which are exogenous and list those in 'interact_exogenous', then identify 'interact_with_exogneous'
   #   
   #   ## something like the below but not quite: 
   #   #interact_exogenous <- exogenous[regexpr(ml, exogenous)<0]
   #  # interact_with_exogenous <- ml[regexpr(ml, exogenous)<0]
   #   
   #   } else{
   #     interatct_exogenous <- NULL
   #     interact_with_exogenous <- NULL
   #     }
   #    
   #   multiLASSO(data                       = data,
   #              out                        = out,
   #              sep                        = sep,
   #              header                     = header,
   #              ar                         = ar,
   #              plot                       = plot,
   #              conv_vars                  = conv_vars,
   #              conv_length                = conv_length,
   #              conv_interval              = conv_interval,
   #              #mult_vars                  = mult_vars,
   #              # mean_center_mult          = FALSE,
   #              # standardize               = FALSE,
   #              groupcutoff                = groupcutoff,
   #              alpha                      = .5,
   #              model_crit                 = lasso_model_crit,
   #              penalties                  = NULL,
   #              test_penalties             = FALSE,
   #              exogenous                  = exogenous,
   #              lag_exogenous              = FALSE,
   #              interact_exogenous         = interact_exogenous, 
   #              interact_with_exogenous    = interact_with_exogenous,
   #              predict_with_interactions  = NULL)
   #   
   # } else {

  dat         <- setup(data                 = data,
                       sep                  = sep,
                       header               = header,
                       out                  = out,
                       plot                 = plot,
                       ar                   = ar,
                       paths                = paths,
                       exogenous            = exogenous,
                       mult_vars            = mult_vars,
                       mean_center_mult     = mean_center_mult,
                       standardize          = standardize,
                       subgroup             = subgroup,
                       ind                  = FALSE,
                       agg                  = FALSE,
                       groupcutoff          = groupcutoff,
                       subcutoff            = subcutoff,
                       conv_vars            = conv_vars, 
                       conv_length          = conv_length, 
                       conv_interval        = conv_interval,
                       lv_model             = lv_model,
                       lv_estimator         = lv_estimator,
                       lv_scores            = lv_scores,
                       lv_miiv_scaling      = lv_miiv_scaling)

  
  #Error Check for Confirm Subgroup Labels
  if(subgroup & !is.null(confirm_subgroup)){
    if(dim(confirm_subgroup)[[2]] != 2){
      stop(paste0("gimme ERROR: confirmatory subgroup dataframe is not a two column dataframe.",
                  " Please ensure that the confirmatory subgroup dataframe consists of a column of filenames and a column of community assignments."))
    }
    if(length(match(confirm_subgroup[,1], (dat$file_order)$names)) != dim(confirm_subgroup)[[1]]){
      stop(paste0("gimme ERROR: confirmatory subgroup dataframe contains mismatched filenames.",
                  " Please ensure that the confirmatory subgroup filenames match the data filenames, sans extensions (Example: sub_1000, not sub_1000.csv)"))
    }
    if(!is.numeric(confirm_subgroup[,2])){
      stop(paste0("gimme ERROR: confirmatory subgroup assignments are non-numeric.",
                  " Please ensure that the confirmatory subgroup assignments are integer valued, beginning from 1. (Example: 1, 2, 3, 4)"))
    }
  }
  
  
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  # Group Search Stage
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  grp <- list(
    "n_group_paths" = 0,
    "n_fixed_paths" = length(dat$fixed_paths),
    "group_paths"   = c()
  )

  if(!hybrid){
    elig_paths   = dat$candidate_paths
  } else{
    elig_paths   = c(dat$candidate_paths, dat$candidate_corr)
  }
  
  grp_hist  <- search.paths(
    base_syntax    = dat$syntax,
    fixed_syntax   = NULL,
    add_syntax     = grp$group_paths,
    n_paths        = grp$n_group_paths,
    data_list      = dat$ts_list,
    elig_paths     = elig_paths, 
    prop_cutoff    = dat$group_cutoff,
    n_subj         = dat$n_subj,
    chisq_cutoff   = qchisq(1-.05/dat$n_subj, 1),
    subgroup_stage = FALSE,
    ms_allow       = ms_allow,
    ms_tol         = ms_tol, 
    hybrid         = hybrid,
    dir_prop_cutoff = dir_prop_cutoff
  )
  
  
  
  #-------------------------------------------------------------#
  # Create group object from grp_hist
  #-------------------------------------------------------------#
  
  grp <- replicate(length(grp_hist[[length(grp_hist)]]), grp, simplify = FALSE)
  
  grp <- lapply(seq_along(grp), function(i){
    grp[[i]]$n_group_paths <- grp_hist[[length(grp_hist)]][[i]]$n_paths
    grp[[i]]$group_paths   <- grp_hist[[length(grp_hist)]][[i]]$add_syntax
    grp[[i]]
  })
   
  
  #-------------------------------------------------------------#
  # Add fields we will need later to grp_hist
  #-------------------------------------------------------------#
  grp_hist <- lapply(seq_along(grp_hist), function(i){
     lapply(seq_along(grp_hist[[i]]), function(j){
       grp_hist[[i]][[j]]$pruned     <- NA
       grp_hist[[i]][[j]]$stage      <- "group"
       grp_hist[[i]][[j]]$grp_sol    <- NA
       grp_hist[[i]][[j]]$sub_sol    <- NA
       grp_hist[[i]][[j]]$ind_sol    <- NA
       if(i == length(grp_hist[[length(grp_hist)]])){
         grp_hist[[i]][[j]]$grp_sol <- j
       }
       grp_hist[[i]][[j]]
     })
  })

  
  
  
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  # Group Prune Stage
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  #-------------------------------------------------------------#
  
  prune <- unlist(lapply(grp,"[","n_group_paths")) != 0

  if (any(prune)){
    
    grp_prune <- lapply(seq_along(grp), function(i){
      
      if(ms_allow){
        writeLines(paste0("solution ", i ,"..."))
      }
      
      prune.paths(
        base_syntax    = dat$syntax,
        fixed_syntax   = NULL,
        add_syntax     = grp[[i]]$group_paths,
        data_list      = dat$ts_list,
        n_paths        = grp[[i]]$n_group_paths,
        n_subj         = dat$n_subj,
        prop_cutoff    = dat$group_cutoff,
        elig_paths     = grp$group_paths,
        subgroup_stage = FALSE
      )
      
    })

  
    
    # #-----------------------------------------------------------#
    # # Add the pruning information to the grp history.
    # Note: the group history does not contain "accurate"
    #       accounting of the model state, as pruned paths
    #       are indicated in the pruned field of the final
    #       entry of a given stage.
    # #-----------------------------------------------------------#

    grp_hist[[length(grp_hist)]] <- lapply(seq_along(grp_hist[[length(grp_hist)]]), function(i){
      pruned <- setdiff(grp[[i]]$group_paths, grp_prune[[i]]$add_syntax)
      if(length(pruned) != 0){
        grp_hist[[length(grp_hist)]][[i]]$pruned <- pruned
      }
      grp_hist[[length(grp_hist)]][[i]]
    })
    
    
    #-----------------------------------------------------------#
    # Now, update the grp object to reflect the pruning.
    #-----------------------------------------------------------#
    
    grp <- lapply(seq_along(grp), function(i){

      grp[[i]][["n_group_paths"]] <- grp_prune[[i]]$n_paths
      grp[[i]][["group_paths"  ]] <- grp_prune[[i]]$add_syntax
      grp[[i]]

    })
    
  } # end pruning 

  
  #-----------------------------------------------------------#
  # Subgroup stage.
  #-----------------------------------------------------------#
  if (subgroup){
  
    sub_res <- lapply(seq_along(grp), function(i){
      res <- subgroupStage(
        dat,
        grp[[i]],
        confirm_subgroup,
        elig_paths,
        sub_feature,
        sub_method,
        ms_tol   = ms_tol,
        ms_allow = FALSE,
        sub_sim_thresh = sub_sim_thresh,
        hybrid, 
        dir_prop_cutoff = dir_prop_cutoff
      )
    })

    
    #-------------------------------------------------------------#
    # Break up the results from subgrouping.
    #-------------------------------------------------------------#
    sub      <- lapply(sub_res, "[[", "sub")
    sub_spec <- lapply(sub_res, "[[", "sub_spec")
    ind      <- lapply(sub_res, "[[", "ind")
    grp_sub  <- lapply(sub_res, "[[", "grp")
    
    
    #-------------------------------------------------------------#
    # Expand the grp_history object to contain subgroup solutions.
    #-------------------------------------------------------------#
    sub_hist  <- list(); cnt<-1
    
    for(i in 1:length(grp)){
      
      sub_pruned_grp_paths <- setdiff(grp[[i]]$group_paths, grp_sub[[i]]$group_paths)
      grp_level_paths      <- grp_hist[[length(grp_hist)]][[i]]$add_syntax
      grp_level_n_paths    <- grp_hist[[length(grp_hist)]][[i]]$n_paths
      
      for(j in 1:sub[[i]]$n_subgroups){
        
        sol_i_subgrp_j_hist            <- grp_hist[[length(grp_hist)]][[i]]
        sol_i_subgrp_j_hist$add_syntax <- c(grp_level_paths, sub_spec[[i]][[j]]$sub_paths)
        sol_i_subgrp_j_hist$n_paths    <- grp_level_n_paths + sub_spec[[i]][[j]]$n_sub_paths
        sol_i_subgrp_j_hist$pruned     <- sub_pruned_grp_paths
        sol_i_subgrp_j_hist$stage      <- "subgroup"
        sol_i_subgrp_j_hist$grp_sol    <- i
        sol_i_subgrp_j_hist$sub_sol    <- j
        sol_i_subgrp_j_hist$ind_sol    <- NA
        
        sub_hist[[cnt]] <- sol_i_subgrp_j_hist
        cnt <- cnt + 1
        
      }
      
    }
    
    hist <- append(grp_hist, grp_hist[length(grp_hist)])
    hist[[length(hist)]] <- sub_hist

    #-------------------------------------------------------------#
    # Update the grp object in case grp paths were pruned
    #-------------------------------------------------------------#
    grp <- grp_sub
    
  } else {

    sub      <- NULL
    sub_spec <- NULL
    ind      <- dat$file_order
    hist     <- grp_hist
    
  }
  
  
  #-------------------------------------------------------------#
  # If this is classic gimme...
  #-------------------------------------------------------------#
  if(!ms_allow){
  
    # 2.19.2019 kmg: ind[1]$ returns NULL for subgroups; changed to ind[[1]] here
    if(subgroup){
      store <- indiv.search(dat, grp[[1]], ind[[1]], hybrid)
    } else {
      store <- indiv.search(dat, grp[[1]], ind[1], hybrid)
    }
    
    if(!is.null(lv_model)){
      dat$lvgimme$miiv_est_table <- fitFinalGimmeModels(
          ts_list_obs = dat$lvgimme$ts_list_obs,
          meas_model  = dat$lvgimme$model_list_dfa,
          lv_model    = lapply(store$syntax, function(x){x[!grepl("0\\*", x)]}),
          miiv.dir    =  if(is.null(dat$out)) NULL else {file.path(dat$out,"miiv")},
          lv_final_estimator = lv_final_estimator
      )
    }
  
    print.gimme(x = sub[[1]],
                y = subgroup,
                z = dat)
  
    # wrap-up and create output
    final <- final.org(dat,
                       grp = grp[[1]],
                       sub = sub[[1]],
                       sub_spec = sub_spec[[1]],
                       diagnos = diagnos,
                       store)
  
    # plot.gimmep convenience functions. 
    # if you change an object name here, 
    # please check plot.gimmep and print.gimmep
    # to ensure compatibility
    
    res <- list(data            = dat$ts_list,
                path_est_mats   = store$betas,
                varnames        = dat$varnames,
                n_vars_total    = dat$n_vars_total,
                n_lagged        = dat$n_lagged,
                n_endog         = dat$n_endog,
                fit             = final$fit,
                path_se_est     = final$param_est,
                plots           = store$plots,
                group_plot_paths= final$samp_plot,
                group_plot_cov  = final$samp_plot_cov,
                sub_plots_paths = final$sub_plots,
                sub_plots_cov   = final$sub_plots_cov,
                subgroup        = subgroup,
                path_counts     = final$sample_counts,
                path_counts_sub = final$sub_counts,
                cov_counts      = final$sample_counts_cov,
                cov_counts_sub  = final$sub_counts_cov,
                vcov            = store$vcov,
                vcovfull        = store$vcovfull,
                psi             = store$psi,
                psi_unstd       = store$psiunstd,
                sim_matrix      = sub[[1]]$sim, 
                syntax          = dat$syntax,
                lvgimme         = dat$lvgimme
                )
    class(res) <- "gimmep"
    
    
    invisible(res)
  
    
  } else {


    #-----------------------------------------------#
    # GIMME-MS: Individual level search             #
    #-----------------------------------------------#
    
    ind_results <- lapply(seq_along(grp), function(j){

      if(ms_allow){
        writeLines(paste0("group solution ", j ,"..."))
      }

      if(subgroup){
        indiv.search.ms(dat, grp[[j]], ind[[j]], ms_tol, ms_allow, j, hybrid = FALSE)
      } else {
        indiv.search.ms(dat, grp[[j]], ind, ms_tol, ms_allow, j, hybrid = FALSE)
      }

    })
    
    ind_hist <- lapply(ind_results, "[[", "ind_history")
    ind_fit  <- lapply(ind_results, "[[", "res")
    
    writeLines("gimme multiple solutions finished running normally")
    
    res <- list(
      dat = dat,
      grp = grp,
      ind = ind,
      sub = sub,
      sub_spec = sub_spec,
      grp_hist = hist, 
      ind_hist = ind_hist,
      ind_fit  = ind_fit
    )
    
    class(res) <- "gimmemsp"
    
    # write gimme ms results to csv file
    res$tables <- gimmems.write(res)
    
    invisible(res)
  }
  
}

print.gimme <- function(x, y, z){
  writeLines("gimme finished running normally")
  if (!is.null(z$out)) writeLines(paste("output is stored in", z$out))
  if (y == TRUE) {
    writeLines(paste("Number of subgroups =", x$n_subgroups))
    writeLines(paste("Modularity =", round(x$modularity, digits = 5)))
  }
}

