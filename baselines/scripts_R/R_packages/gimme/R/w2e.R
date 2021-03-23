#' Create edge list from weight matrix.
#' @param x The coefficient matrix from an individual
#' @return A list of all non-zero edges to feed to qgraph
#' @keywords internal
w2e <- function(x) cbind(which(x != 0, arr.ind = TRUE), x[x != 0])