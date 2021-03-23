#' @name indSEM
#' @aliases indSEM
#' @title Individual-level structural equation model search.
#' @description This function identifies structural equation models for each
#' individual. It does not utilize any shared information from the sample.
#' @usage
#' indSEM(data   = NULL,
#'        out    = NULL,
#'        sep    = NULL,
#'        header = NULL,
#'        ar     = TRUE,
#'        plot   = TRUE,
#'        paths  = NULL,
#'        exogenous        = NULL, 
#'        conv_vars        = NULL,
#'        conv_length      = 16, 
#'        conv_interval    = 1,
#'        mult_vars        = NULL,
#'        mean_center_mult = FALSE,
#'        standardize      = FALSE,
#'        hybrid = FALSE)
#' @param data The path to the directory where the data files are located, 
#' or the name of the list containing each individual's time series. Each file 
#' or matrix must contain one matrix for each individual containing a 
#' T (time) by p (number of variables) matrix where the columns represent 
#' variables and the rows represent time. 
#' @param out The path to the directory where the results will be stored 
#' (optional). If specified, a copy of output files will be replaced in 
#' directory. If directory at specified path does not exist, it will be created.
#' @param sep The spacing of the data files. "" indicates space-delimited, 
#' "/t" indicates tab-delimited, "," indicates comma delimited. Only necessary
#' to specify if reading data in from physical directory.
#' @param header Logical. Indicate TRUE for data files with a header. 
#' Only necessary to specify if reading data in from physical directory.
#' @param ar Logical. If TRUE, begins search for individual models with 
#' autoregressive (AR) paths open. Defaults to TRUE.
#' @param plot Logical. If TRUE, graphs depicting relations among variables of 
#' interest will automatically be created. Defaults to TRUE. For individual-
#' level plots, red paths represent positive weights and blue paths represent
#' negative weights.
#' @param paths lavaan-style syntax containing paths with which to begin model
#' estimation. That is, Y~X indicates that Y is regressed on X, or X 
#' predicts Y. If no header is used, then variables should be referred to with 
#' V followed (with no separation) by the column number. If a header is used, 
#' variables should be referred to using variable names. To reference lag 
#' variables, "lag" should be added to the end of the variable name with no 
#' separation. Defaults to NULL.
#' @param exogenous Vector of variable names to be treated as exogenous.  
#' That is, exogenous variable X can predict Y  but cannot be predicted by Y.  
#' If no header is used, then variables should be referred to with V followed 
#' (with no separation) by the column number. If a header is used, variables 
#' should be referred to using variable names.  Defaults to NULL.
#' @param conv_vars Vector of variable names to be convolved via smoothed Finite Impulse 
#' Response (sFIR). Defaults to NULL.
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
#' variable name should be followed by "lag" with no separation (e.g. V1*V2lag). Note that if
#' multiplied variables are desired, at least one variable in the dataset must be specified as exogenous.
#' Defaults to NULL.
#' @param mean_center_mult Logical. If TRUE, the variables indicated in mult_vars will be mean-centered
#' before being multiplied together. Defaults to FALSE. 
#' @param standardize Logical. If TRUE, all variables will be standardized to have a mean of zero and a
#' standard deviation of one. Defaults to FALSE. 
#' @param hybrid Logical. If TRUE, enables hybrid-VAR models where both directed contemporaneous paths and contemporaneous 	
#' covariances among residuals are candidate relations in the search space. Defaults to FALSE.
#' @details
#'  In main output directory:
#'  \itemize{
#'  \item{\strong{indivPathEstimates}} {Contains estimate, standard error, 
#'  p-value, and z-value for each path for each individual}
#'  \item{\strong{summaryFit}} {Contains model fit information for individual-
#'  level models. }
#'  \item{\strong{summaryPathCountMatrix}} Contains counts of total number of 
#'  paths, both contemporaneous and lagged, estimated for the sample. The row
#'   variable is the outcome and the column variable is the predictor variable.
#'  \item{\strong{summaryPathCounts}} {Contains summary count information for 
#'  paths identified at the individual-level.}
#'  \item{\strong{summaryPathsPlot}} Contains counts of total number of paths,
#'   both contemporaneous and lagged, estimated for the sample. The row variable
#'    is the outcome and the column variable is the predictor variable.
#'  }
#'  In individual output directory (where \strong{\emph{id}} represents the
#'   original file name for each individual):
#'  \itemize{
#'  \item{\strong{\emph{id}Betas}} Contains individual-level estimates of each 
#'  path for each individual.
#'  \item{\strong{\emph{id}StdErrors}} Contains individual-level standard errors 
#'  for each path for each individual.
#'  \item{\strong{\emph{id}Plot}} Contains individual-level plots. Red paths 
#'  represent positive weights and blue paths represent negative weights.
#' }
#' @author Stephanie Lane
#' @examples
#'  \dontrun{
#' fit <- indSEM(data   = "C:/data100",
#'               out    = "C:/data100_indSEM_out",
#'               sep    = ",",
#'               header = FALSE)
#' print(fit, file = "group1.1", estimates = TRUE)
#' plot(fit, file = "group1.1")
#'  }
#'@keywords indSEM
#'@export

indSEM <- function(data   = NULL,
                   out    = NULL,
                   sep    = NULL,
                   header = NULL,
                   ar     = TRUE,
                   plot   = TRUE,
                   paths  = NULL,
                   exogenous = NULL, 
                   conv_vars      = NULL,
                   conv_length    = 16, 
                   conv_interval = 1, 
                   mult_vars      = NULL,
                   mean_center_mult = FALSE,
                   standardize    = FALSE,
                   hybrid = FALSE){
  
  dat  <- setup(data        = data,
                sep         = sep,
                header      = header,
                out         = out,
                plot        = plot,
                ar          = ar,
                paths       = paths,
                exogenous   = exogenous,
                mult_vars   = mult_vars,
                mean_center_mult  = mean_center_mult,  
                standardize = standardize,
                conv_vars = conv_vars, 
                conv_length = conv_length, 
                conv_interval = conv_interval,
                groupcutoff = NULL,
                subcutoff   = NULL,
                subgroup    = FALSE,
                ind         = TRUE,
                agg         = FALSE)

  store <- indiv.search(dat, 
                        grp = NULL, 
                        ind = dat$file_order,
                        hybrid)

  final <- final.org(dat, 
                     grp      = NULL, 
                     sub      = NULL, 
                     sub_spec = NULL, 
                     diagnos = FALSE,
                     store = store)

  print.gimme.indSEM(z = dat)
  
  res <- list(path_est_mats = store$betas,
              varnames = dat$varnames,
              n_rois = dat$n_rois,
              fit = final$fit,
              path_se_est = final$param_est,
              plots = store$plots,
              group_plot = final$samp_plot,
              group_plot_cov = final$samp_plot_cov,
              path_counts = final$sample_counts,
              cov_counts      = final$sample_counts_cov,
              vcov = store$vcov,
              vcovfull = store$vcovfull)
  
  class(res) <-  "indSEMp"
  
  invisible(res)
}

print.gimme.indSEM <- function(z){
  writeLines("indSEM finished running normally")
  if (!is.null(z$out)) writeLines(paste("output is stored in", z$out))
}
