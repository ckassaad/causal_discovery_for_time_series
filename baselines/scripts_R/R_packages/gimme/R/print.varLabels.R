#' Print method for varLabels object.
#' @keywords internal
print.varLabels  <- function(x){

  cat(
    " Original Variable Names \n  ", paste0(x$orig, collapse = ", "), "\n\n",
    "User-Specified Exogenous Variable Names \n  ", paste0(x$uexo, collapse = ", "), "\n\n",
    "Lagged Variable Names \n  ", paste0(x$lagg, collapse = ", "), "\n\n",
    "Convolved Variable Names \n  ", paste0(x$conv, collapse = ", "), "\n\n",
    "Multiplied Variable Names \n  ", paste0(x$mult, collapse = ", "), "\n\n",
    "Exogenous Variable Names \n  ", paste0(x$exog, collapse = ", "), "\n\n",
    "Endogenous Variable Names \n  ", paste0(x$endo, collapse = ", "), "\n\n",
    "Categorical Variable Names \n  ", paste0(x$catg, collapse = ", "), "\n\n",
    "Standardized Variable Names \n  ", paste0(x$stnd, collapse = ", "), "\n\n",
    "All Variable Names (column order) \n  ", paste0(x$coln, collapse = ", "), "\n\n"
  )
  
}


