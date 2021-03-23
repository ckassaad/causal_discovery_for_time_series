#' Attempt to fit lavaan model.
#' @param syntax A character vector containing syntax.
#' @param data_file A data frame containing individual data set.
#' @return If successful, returns fitted lavaan object. If not successful,
#' catches and returns error message.
#' @keywords internal 
fit.model <- function (syntax, data_file) {
  
  fit <- tryCatch(lavaan::lavaan(syntax,
                                 data            = data_file,
                                 model.type      = "sem",
                                 missing         = "fiml",
                                 estimator       = "ml",
                                 int.ov.free     = FALSE,
                                 int.lv.free     = TRUE,
                                 auto.fix.first  = TRUE,
                                 auto.var        = TRUE,
                                 auto.cov.lv.x   = TRUE,
                                 auto.th         = TRUE,
                                 auto.delta      = TRUE,
                                 auto.cov.y      = FALSE,
                                 auto.fix.single = TRUE,
                                 warn            = FALSE),
                  error = function(e) e)
  return(fit)
}