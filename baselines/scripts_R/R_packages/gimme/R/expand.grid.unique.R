#' Provides unique combinations of two vectors.
#' @param x A character vector containing variable names.
#' @param y A character vector containing variable names.
#' @param incl.eq Logical. TRUE means that combinations are kept where
#' a variable appears twice.
#' @return The unique combinations of the variable names. Used in syntax
#' creation.
#' @keywords internal 
expand.grid.unique <- function(x, y, incl.eq = TRUE){
  g <- function(i){
    z <- setdiff(y, x[seq_len(i - incl.eq)])
    if(length(z)) cbind(x[i], z, deparse.level = 0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}