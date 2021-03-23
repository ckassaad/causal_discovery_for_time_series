#' Create structure of group search solutions.
#' @keywords internal 
create.tree <- function(history, subgroup, individual = FALSE, all.ind = FALSE){
  
  ###
  # history <- ind_hist_k; grp_hist_i
  # individual <- FALSE
  ##     
       
  
  #-------------------------------------------------------#
  # Internal functions for tree building
  #-------------------------------------------------------#
  rbind.NA.fill <- function (...) {
    dargs <- list(...)
    if (!all(vapply(dargs, is.vector, TRUE)))
      stop("all inputs must be vectors")
    if (!all(vapply(dargs, function(x) !is.null(names(x)), TRUE)))
      stop("all input vectors must be named.")
    all.names <- unique(names(unlist(dargs)))
    out <- do.call(rbind, lapply(dargs, `[`, all.names))
    colnames(out) <- all.names
    out
  }

  paste4 <- function(x, sep = ", ") {
    x <- gsub("^\\s+|\\s+$", "", x)
    ret <- paste(x[!is.na(x) & !(x %in% "")], collapse = sep)
    is.na(ret) <- ret == ""
    return(ret)
  }
  
  #-------------------------------------------------------#
  # Formatting for the additional columns
  #-------------------------------------------------------#
  format_character_vars <- function(x){
    if (!is.na(x)) {
      res <- as.character(x) 
    } else {
      res <- "" 
    } 
    return(res)
  }
      
  format_number_vars <- function(x){
    if (!is.na(x)) {
      res <- round(x,3)
    } else {
      res <- "" 
    } 
    return(res)
  }
  
  format_sub_sol <- function(x){
    if (is.na(x)) {
      res <- ""
    } else {
      if(length(x)>1){
        res <- paste0(x,collapse=",")
      } else {
        res <- x; #as.character(x)
      }
    
    } 
    return(res)
  }

  
  
  #-------------------------------------------------------#
  # For the group model
  #-------------------------------------------------------#
  if(!individual){
    
    if(!subgroup){
    
      hist_list <- lapply(seq_along(history), function(i){ 
        history[[i]]$add_syntax
      })
      
      new_list <- list()
      
      for(i in 1:length(hist_list)){ 
        tmp_list <- hist_list[[i]]
        names(tmp_list) <- paste0("V", 1:length(tmp_list))
        new_list <- append(new_list, list(tmp_list))
      }
  
      hist_list <- new_list
      
      #-------------------------------------------------------#
      # Create variables to go beside tree plot
      #-------------------------------------------------------#
      stage  <- unlist(lapply(history, function(x){x$stage}))
      grp_sol <- unlist(lapply(history, function(x){x$grp_sol}))
      pruned  <- unlist(lapply(seq_along(history), function(i){ 
        pruned <- history[[i]]$pruned
        if(any(is.na(pruned))){ 
          NA 
        } else {
          paste0(pruned,collapse=", ")
        }
      }))  
      #-------------------------------------------------------#
      # Create tree
      #-------------------------------------------------------#
      df <- as.data.frame(do.call(rbind.NA.fill, hist_list))
      V0 <- rep("NULL MODEL", nrow(df))
      df <- cbind(V0, df)
      df$pathString <- apply(df, 1, function(x) paste4(x, sep = "/"))
      df$grp_sol <- grp_sol
      df$stage   <- stage
      df$pruned <- pruned
          
      #-------------------------------------------------------#
      # Format columns
      #-------------------------------------------------------#
      dt <- data.tree::as.Node(df)
      data.tree::SetFormat(dt, "stage", format_character_vars)
      data.tree::SetFormat(dt, "grp_sol", format_character_vars)
      data.tree::SetFormat(dt, "pruned", format_character_vars)
   
    #-------------------------------------------------------#
    # For the group and subgroup model
    #-------------------------------------------------------#
    } else {
      
      hist_list <- lapply(seq_along(history), function(i){ 
        lapply(seq_along(history[[i]]), function(j){ 
          history[[i]][[j]]$add_syntax
        })
      })
      
      new_list <- list()
      for(i in 1:length(hist_list)){ 
        for(j in 1:length(hist_list[[i]])){
          tmp_list <- hist_list[[i]][[j]]
          names(tmp_list) <- paste0("V", 1:length(tmp_list))
          new_list <- append(new_list, list(tmp_list))
        }
      }
  
      hist_list <- new_list
      
      #-------------------------------------------------------#
      # Create variables to go beside tree plot
      #-------------------------------------------------------#
      stage  <- unlist(lapply(history, function(x){
        lapply(x, "[[", "stage")
      }))
      grp_sol <- unlist(lapply(history, function(x){
        lapply(x, "[[", "grp_sol")
      }))
      sub_sol <- unlist(lapply(history, function(x){
        lapply(x, "[[", "sub_sol")
      }))
      modularity <- unlist(lapply(history, function(x){
        lapply(x, "[[", "modularity")
      }))
      pruned  <- unlist(lapply(seq_along(history), function(i){ 
        lapply(seq_along(history[[i]]), function(j){
          #print(history[[i]][[j]]$pruned)
          if(length( history[[i]][[j]]$pruned)>0){
            if(any(is.na( history[[i]][[j]]$pruned))){
              NA
            } else {
              paste0( history[[i]][[j]]$pruned,collapse=", ")
            }
          } else {
            NA
          }
        })
      }))  
      #-------------------------------------------------------#
      # Create tree
      #-------------------------------------------------------#
      df <- as.data.frame(do.call(rbind.NA.fill, hist_list))
      V0 <- rep("NULL MODEL", nrow(df))
      df <- cbind(V0, df)
      df$pathString <- apply(df, 1, function(x) paste4(x, sep = "/"))
      df$grp_sol <- grp_sol
      df$sub_sol <- sub_sol
      df$modulatiry <- modularity
      df$stage   <- stage
      df$pruned <- pruned
      
      foo <- function(x){as.vector(unique(x))}
      
      sub_sol_df <- aggregate(
        sub_sol ~ pathString, 
        FUN=foo,
        data = df,
        na.action=NULL
      )
      
      grp_sol_df <- aggregate(
        grp_sol ~ pathString, 
        FUN=foo,
        data = df,
        na.action=NULL
      )
      
      df.tmp <- merge(sub_sol_df,grp_sol_df,by="pathString")
      
      mod_df <- aggregate(
        modularity ~ pathString, 
        FUN=foo,
        data = df,
        na.action=NULL
      )
      
      df.tmp <- merge(df.tmp,mod_df,by="pathString")
      
      pruned_df <- aggregate(
        pruned ~ pathString, 
        FUN=foo,
        data = df,
        na.action=NULL
      )
      
      df.tmp <- merge(df.tmp,pruned_df,by="pathString")
      
      stage_df <- aggregate(
        stage ~ pathString, 
        FUN=foo,
        data = df,
        na.action=NULL
      )
      
      df.tmp <- merge(df.tmp,stage_df,by="pathString")
  
      #-------------------------------------------------------#
      # Format columns
      #-------------------------------------------------------#
      dt <- data.tree::as.Node(df.tmp)
      data.tree::SetFormat(dt, "stage", format_character_vars)
      data.tree::SetFormat(dt, "grp_sol", format_sub_sol)
      data.tree::SetFormat(dt, "sub_sol", format_sub_sol)
      data.tree::SetFormat(dt, "modularity", format_number_vars)
      data.tree::SetFormat(dt, "pruned", format_character_vars)
      
      #print(dt,"stage","grp_sol","sub_sol","modularity","pruned")
      
    }
   
  #-------------------------------------------------------#
  # Individual Models start
  #-------------------------------------------------------# 
  } else if (individual) {
    
    
    if(!all.ind){
    
      #-------------------------------------------------------#
      # Start building histories
      #-------------------------------------------------------#
      ind_hist <- history
        
      # ind_hist <- ind_hist_k
          
      hist_list <- lapply(seq_along(ind_hist), function(i){ 
        lapply(seq_along(ind_hist[[i]]), function(j){
          lapply(ind_hist[[i]][[j]], "[[", "add_syntax")
        })
      })
        
      new_list <- list()
        
      for(i in 1:length(hist_list)){ 
        for(j in 1:length(hist_list[[i]])){
          for(k in 1:length(hist_list[[i]][[j]])){
            if(!is.null(hist_list[[i]][[j]][[k]])){
              tmp_list <- hist_list[[i]][[j]][[k]]
              names(tmp_list) <- paste0("V", 1:length(tmp_list))
              new_list <- append(new_list, list(tmp_list))
            } else {
              temp_vec <- c(NA)
              names(temp_vec) <- "V1"
              new_list <- append(new_list, list(temp_vec))
            }
          }
        }
      }
          
      hist_list <- new_list
      
      #-------------------------------------------------------#
      # Create variables to go beside tree plot
      #-------------------------------------------------------#
      stage  <- unlist(lapply(ind_hist, function(x){
        lapply(x, function(y){lapply(y, "[[", "stage")})
      }))
          
      grp_sol  <- unlist(lapply(seq_along(ind_hist), function(i){ 
        lapply(seq_along(ind_hist[[i]]), function(j){
          lapply(ind_hist[[i]][[j]], "[[", "group_sol")
        })
      }))
        
      ind_sol  <- unlist(lapply(seq_along(ind_hist), function(i){ 
        lapply(seq_along(ind_hist[[i]]), function(j){
          lapply(ind_hist[[i]][[j]], "[[", "ind_sol")
        })
      }))
        
      pruned  <- unlist(lapply(seq_along(ind_hist), function(i){ 
        lapply(seq_along(ind_hist[[i]]), function(j){
          lapply(seq_along(ind_hist[[i]][[j]]),  function(l){
            pruned <- ind_hist[[i]][[j]][[l]]$pruned
            if(any(is.na(pruned))){ 
              NA 
            } else {
              paste0(pruned,collapse=", ")
            }
          })
        })
      }))
        
      fit <- list(); cnt <- 1
      for(i in 1:length(ind_hist)){
        for(j in 1:length(ind_hist[[i]])){
          for(l in 1:length(ind_hist[[i]][[j]])){
            if(is.na(ind_hist[[i]][[j]][[l]]$fit)[1]){
              fit[[cnt]] <- c("NA" = NA)
            } else {
              fit[[cnt]] <- unclass(ind_hist[[i]][[j]][[l]]$fit)
            }
            cnt <- cnt + 1
          }
        }
      }
    
      fit <- as.data.frame(do.call(rbind.NA.fill, fit))[,-1]
      
      #-------------------------------------------------------#
      # Create tree
      #-------------------------------------------------------#
      df <- as.data.frame(do.call(rbind.NA.fill, hist_list))
      V0 <- rep("NULL MODEL", nrow(df))
      df <- cbind(V0, df)
      df$pathString <- apply(df, 1, function(x) paste4(x, sep = "/"))
        
      
      #-------------------------------------------------------#
      # Individual level plot without subgrouping
      #-------------------------------------------------------#
      if(!subgroup){
        
        df$grp_sol <- grp_sol
        df$stage <- stage
        df$pruned <- pruned
        df$ind_sol <- ind_sol
        fit[is.na(fit)] <- ""
        df <- cbind(df, fit)
          
        #-------------------------------------------------------#
        # Format columns
        #-------------------------------------------------------#
        dt <- data.tree::as.Node(df)
        data.tree::SetFormat(dt, "stage", format_character_vars)
        data.tree::SetFormat(dt, "grp_sol", format_character_vars)
        data.tree::SetFormat(dt, "pruned", format_character_vars)
        data.tree::SetFormat(dt, "ind_sol", format_character_vars)
        
        
      #-------------------------------------------------------#
      # Individual level plot without subgrouping
      #-------------------------------------------------------#   
      } else if(subgroup) {
        
        sub_sol  <- unlist(lapply(seq_along(ind_hist), function(i){ 
          lapply(seq_along(ind_hist[[i]]), function(j){
            lapply(ind_hist[[i]][[j]], "[[", "sub_sol")
          })
        }))
        
    
        df$grp_sol <- grp_sol
        df$sub_sol <- sub_sol
        df$stage   <- stage
        df$pruned  <- pruned
        df$ind_sol <- ind_sol
        
        if(individual){
          
          fit[is.na(fit)] <- ""
          df <- cbind(df, fit)
          
        } else {
          
          df <- df[df$stage == "group" | df$stage == "subgroup",]
          
        }
        
        #-------------------------------------------------------#
        # Format columns
        #-------------------------------------------------------#
        dt <- data.tree::as.Node(df)
        data.tree::SetFormat(dt, "stage", format_character_vars)
        data.tree::SetFormat(dt, "grp_sol", format_character_vars)
        data.tree::SetFormat(dt, "sub_sol", format_character_vars)
        data.tree::SetFormat(dt, "pruned", format_character_vars)
        data.tree::SetFormat(dt, "ind_sol", format_character_vars)
      
        
      }
    
    } else if (all.ind){
      
      # #-------------------------------------------------------#
      # # Start building histories
      # #-------------------------------------------------------#
      # ind_hist <- history
      #   
      # # ind_hist <- ind_hist_k
      #     
      # hist_list <- lapply(seq_along(ind_hist), function(i){ 
      #   lapply(seq_along(ind_hist[[i]]), function(j){
      #     lapply(seq_along(ind_hist[[i]][[j]]), function(l){
      #       lapply(ind_hist[[i]][[j]][[l]], "[[", "add_syntax")
      #     })
      #   })
      # })
      #   
      # new_list <- list()
      #   
      # for(i in 1:length(hist_list)){ 
      #   for(j in 1:length(hist_list[[i]])){
      #     for(k in 1:length(hist_list[[i]][[j]])){
      #       for(l in 1:length(hist_list[[i]][[j]][[k]])){
      #         tmp_list <- hist_list[[i]][[j]][[k]][[l]]
      #         names(tmp_list) <- paste0("V", 1:length(tmp_list))
      #         new_list <- append(new_list, list(tmp_list))
      #       }
      #     }
      #   }
      # }
      #     
      # hist_list <- new_list
      # 
      # #-------------------------------------------------------#
      # # Create variables to go beside tree plot
      # #-------------------------------------------------------#
      # 
      # stage  <- unlist(lapply(seq_along(ind_hist), function(i){ 
      #   lapply(seq_along(ind_hist[[i]]), function(j){
      #     lapply(seq_along(ind_hist[[i]][[j]]), function(l){
      #       lapply(ind_hist[[i]][[j]][[l]], "[[", "stage")
      #     })
      #   })
      # }))
      #     
      # grp_sol  <- unlist(lapply(seq_along(ind_hist), function(i){ 
      #   lapply(seq_along(ind_hist[[i]]), function(j){
      #     lapply(seq_along(ind_hist[[i]][[j]]), function(l){
      #       lapply(ind_hist[[i]][[j]][[l]], "[[", "group_sol")
      #     })
      #   })
      # }))
      #   
      # ind_sol  <- unlist(lapply(seq_along(ind_hist), function(i){ 
      #   lapply(seq_along(ind_hist[[i]]), function(j){
      #     lapply(seq_along(ind_hist[[i]][[j]]), function(l){
      #       lapply(ind_hist[[i]][[j]][[l]], "[[", "ind_sol")
      #     })
      #   })
      # }))
      # 
      # subj_id  <- unlist(lapply(seq_along(ind_hist), function(i){ 
      #   lapply(seq_along(ind_hist[[i]]), function(j){
      #     lapply(seq_along(ind_hist[[i]][[j]]), function(l){
      #       lapply(ind_hist[[i]][[j]][[l]], "[[", "subj_id")
      #     })
      #   })
      # }))
      #   
      # pruned  <- unlist(lapply(seq_along(ind_hist), function(i){ 
      #   lapply(seq_along(ind_hist[[i]]), function(j){
      #     lapply(seq_along(ind_hist[[i]][[j]]), function(k){
      #       lapply(seq_along(ind_hist[[i]][[j]][[k]]),  function(l){
      #         pruned <- ind_hist[[i]][[j]][[k]][[l]]$pruned
      #           if(length(pruned)>0){
      #             if(any(is.na(pruned))){
      #               NA
      #             } else {
      #               paste0(pruned,collapse=", ")
      #             }
      #           } else {
      #             NA
      #           }
      #       })
      #     })
      #   })
      # }))
      #   
      # #-------------------------------------------------------#
      # # Create tree
      # #-------------------------------------------------------#
      # df <- as.data.frame(do.call(rbind.NA.fill, hist_list))
      # V0 <- rep("NULL MODEL", nrow(df))
      # df <- cbind(V0, df)
      # df$pathString <- apply(df, 1, function(x) paste4(x, sep = "/"))
      # df$grp_sol <- grp_sol
      # df$ind_sol <- ind_sol
      # df$stage   <- stage
      # df$pruned  <- pruned
      # df$subj_id  <- subj_id
      # 
      # foo <- function(x){as.vector(unique(x))}
      # 
      # ind_sol_df <- aggregate(
      #   ind_sol ~ pathString, 
      #   FUN=foo,
      #   data = df,
      #   na.action=NULL
      # )
      # 
      # grp_sol_df <- aggregate(
      #   grp_sol ~ pathString, 
      #   FUN=foo,
      #   data = df,
      #   na.action=NULL
      # )
      # 
      # df.tmp <- merge(ind_sol_df,grp_sol_df,by="pathString")
      # 
      # subj_df <- aggregate(
      #   subj_id ~ pathString, 
      #   FUN=foo,
      #   data = df,
      #   na.action=NULL
      # )
      # 
      # df.tmp <- merge(df.tmp,subj_df,by="pathString")
      # 
      # pruned_df <- aggregate(
      #   pruned ~ pathString, 
      #   FUN=foo,
      #   data = df,
      #   na.action=NULL
      # )
      # 
      # df.tmp <- merge(df.tmp,pruned_df,by="pathString")
      # 
      # stage_df <- aggregate(
      #   stage ~ pathString, 
      #   FUN=foo,
      #   data = df,
      #   na.action=NULL
      # )
      # 
      # df.tmp <- merge(df.tmp,stage_df,by="pathString")
      # 
      # #-------------------------------------------------------#
      # # Format columns
      # #-------------------------------------------------------#
      # dt <- data.tree::as.Node(df.tmp)
      # data.tree::SetFormat(dt, "stage", format_character_vars)
      # data.tree::SetFormat(dt, "grp_sol", format_sub_sol)
      # data.tree::SetFormat(dt, "ind_sol", format_sub_sol)
      # data.tree::SetFormat(dt, "subj_id", format_sub_sol)
      # data.tree::SetFormat(dt, "pruned", format_character_vars)
      
      #print(dt,"stage","grp_sol","ind_sol","pruned","subj_id")
      
    } # end all.ind
    
    
  } # end individual
  
  return(dt)

}
  