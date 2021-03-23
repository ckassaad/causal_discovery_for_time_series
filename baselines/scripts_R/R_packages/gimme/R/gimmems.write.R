#' Write MS-GIMME results to data.frame.
#' @keywords internal 
gimmems.write <- function(x){
    
  out.dir <- x$dat$ctrlOpts$out
  
  #------------------------------------------#
  # Summary fit
  #------------------------------------------#
  
  summaryFit <- as.data.frame(do.call("rbind", 
    
    # grpsol_i                       
    lapply(seq_along(x$ind_fit), function(i){                
      
      as.data.frame(do.call("rbind",
                            
        # grpsol_i_ind_j
        lapply(seq_along(x$ind_fit[[i]]), function(j){         
          
          as.data.frame(do.call("rbind",
                                        
            # grpsol_i_ind_j_indsol_k
            lapply(seq_along(x$ind_fit[[i]][[j]]), function(k){  
            
              grp_sol_i_ind_j_sol_k <-  x$ind_fit[[i]][[j]][[k]]
              
              # we want to add some identifying information
              id_cols <- data.frame(
                "subj" = grp_sol_i_ind_j_sol_k$subj,
                "grp_sol" = i,
                "ind_sol" = k
              )
              
              # check fit info from lavaan
              
              fit_info_lav <- grp_sol_i_ind_j_sol_k$fits
              
              # are the elements NA
              
              if(all(is.na(fit_info_lav))){
                
                # zf: designed this way so gimme will throw
                #     an error here if the fit object changes.
                #     this patch is needed because
                #     gimme passes unnamed vector of NAs 
                #     for fit info when model does not converge.
                names(fit_info_lav) <- c(
                  "chisq", "df", "npar", "pvalue", "rmsea", 
                  "srmr", "nnfi", "cfi", "bic", "aic", "logl"
                )
                
              }
                            
              df <- cbind(id_cols,  t(data.frame(fit_info_lav)))
              
              df$status <- grp_sol_i_ind_j_sol_k$status
              
              if( df$status == "converged normally"){
                
                diag.psi <- diag(grp_sol_i_ind_j_sol_k$psi)
                
                improper.var  <- any(diag.psi > 1 | diag.psi < 0)
                
                df$checkPsi <- improper.var 
                
              } else {
                
                df$checkPsi <- NA
                
              }
              
                
              df
            
            })
            
          ))
          
        })
      
      ))
      
    })
    
  ))
  
  if(!is.null(out.dir)){
    # now let's write this df to the ms directory
    write.csv(summaryFit, file.path(out.dir, "summaryFit.csv"), row.names=FALSE)
  }
  
  
  indivPathEstimates <- as.data.frame(do.call("rbind", 
    
    # grpsol_i                       
    lapply(seq_along(x$ind_fit), function(i){                
      
      as.data.frame(do.call("rbind",
                            
        # grpsol_i_ind_j
        lapply(seq_along(x$ind_fit[[i]]), function(j){         
          
          as.data.frame(do.call("rbind",
                                        
            # grpsol_i_ind_j_indsol_k
            lapply(seq_along(x$ind_fit[[i]][[j]]), function(k){  
            
              grp_sol_i_ind_j_sol_k <-  x$ind_fit[[i]][[j]][[k]]
              
              df <- grp_sol_i_ind_j_sol_k$coefs
              
              # we can clean up this df a little. for example,
              # we don't want to include any impossible paths.
              # such as those where lag is on the LHS. 
              
              if(nrow(df) == 1 & all(is.na(df[,"lhs"])) & all(is.na(df[,"op"]))){
                
                # do nothing
                
              } else {
                
                 df <- df[!grepl("lag", df[,"lhs"]) & grepl("~", df[,"op"]),]
                
              }
              
              # we want to add some identifying information
              id_cols <- data.frame(
                "subj" = grp_sol_i_ind_j_sol_k$subj,
                "grp_sol" = i,
                "ind_sol" = k
              )
              
              df <- cbind(id_cols, df)
              
              
              df$op <- df$ci.lower <- df$ci.upper <- NULL
              
              
              colnames(df) <-c("subj", "grp_sol", "ind_sol", 
                               "dv", "iv", "beta", "se", 
                               "z", "pval")
              
              df
            
            })
            
          ))
          
        })
      
      ))
      
    })
    
  ))
  
  if(!is.null(out.dir)){
    # now let's write this df to the ms directory
    write.csv(indivPathEstimates, file.path(out.dir, "indivPathEstimates.csv"), row.names=FALSE)
  }
  
  return(list("indivPathEstimates" = indivPathEstimates, "summaryFit" = summaryFit))

}