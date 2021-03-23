#' Determines subgroups.
#' @param base_syntax A character vector containing syntax that never changes.
#' @param data_list A list of all datasets.
#' @param n_subj The number of subjects in the sample.
#' @param file_order A data frame containing the order of the files and 
#' the names of the files. Used to merge in subgroup assignment and preserve
#' order. 
#' @param chisq_cutoff Cutoff used in order for MI to be considered significant.
#' @param elig_paths A character vector containing eligible paths that
#' gimme is allowed to add to the model. Ensures only EPCs from allowable paths
#' are considered in the creation of the similarity matrix.
#' @param confirm_subgroup A dataframe with the first column a string vector of data file names
#' without extensions and the second vector a integer vector of subgroup labels.
#' @return Returns sub object containing similarity matrix, the number of
#' subgroups, the modularity associated with the subgroup memberships, 
#' and a data frame containing the file names and subgroup memberships.
#' @keywords internal 
determine.subgroups <- function(data_list,
                                base_syntax, 
                                n_subj,
                                chisq_cutoff,
                                file_order,
                                elig_paths,
                                confirm_subgroup, 
                                out_path = NULL,
                                sub_feature,
                                sub_method,
                                sub_sim_thresh,
                                hybrid,
                                dir_prop_cutoff){
  #######################
  # base_syntax  = c(dat$syntax, grp[[i]]$group_paths)
  # data_list    = dat$ts_list
  # n_subj       = dat$n_subj
  # chisq_cutoff = dat$chisq_cutoff_mi_epc
  # file_order   = dat$file_order
  # elig_paths   = dat$candidate_paths
  # confirm_subgroup = confirm_subgroup
  # out_path     = dat$out
  # sub_feature  = sub_feature
  #######################
  writeLines(paste0("subgroup-level search"))
  sub_membership  = NULL # appease CRAN check
  
  sub     <- list()
  z_list  <- list()
  mi_list <- list()
  converge <- matrix(,n_subj,1)

  
  fit <- lapply(seq_along(data_list), function(i){fit.model(
    syntax= base_syntax,
    data_file = data_list[[i]])
  })
  for (k in 1:n_subj){
    # writeLines(paste0("subgroup search, subject ", k, " (",names(data_list)[k],")"))
    # fit          <- fit.model(syntax    = base_syntax,
    #                           data_file = data_list[[k]])
    z_list[[k]]  <- return.zs(fit[[k]])
    mi_list[[k]] <- return.mis(fit[[k]], elig_paths)
    converge[k]  <- lavInspect(fit[[k]], "converged")
  }
  
  names(z_list) <- names(mi_list) <- names(data_list)
  
  # drop individuals who did not converge
  # kmg: this is an imperfect approach because the z_list is NA if no paths 
  # were added. Caused errors, commenting out, added "converge" to return.zs output. 
  # drop    <- unique(c(which(is.na(z_list)), which(is.na(mi_list))))
   drop <- which(converge==FALSE)
   
  if (length(drop) != 0){
    mi_list <- mi_list[-drop]
    z_list  <- z_list[-drop]
  }
  
  mi_list_temp <- lapply(mi_list, 
                         function(x){x$param <- paste0(x$lhs, x$op, x$rhs)
                         x$sig <- x$dir <- ifelse(x$mi > chisq_cutoff, 1, 0)
                         return(x)})
  
  mi_list <- lapply(mi_list_temp, 
                    function(x){subset(x, x$param %in% elig_paths)})
  
  # consider which direction is better when adding contemporaneous similarity; not used in hybrid for now
  if(!hybrid && dir_prop_cutoff>0){
  for (p in 1:length(mi_list)) {
    for (r in 1:length(mi_list[[p]][,1])) {
      if (mi_list[[p]]$dir[r] == 1) {
        grab_opp_mi  <- mi_list[[p]][which(mi_list[[p]]$lhs == mi_list[[p]][r,]$rhs),]
        opp_mi_value <- grab_opp_mi[which(grab_opp_mi$rhs == mi_list[[p]][r,]$lhs),]$mi
        if (length(opp_mi_value)>0) {
          if(opp_mi_value>mi_list[[p]][r,]$mi) 
            mi_list[[p]][r,]$dir <- 0
        }
      }
      opp_mi_value <- 0
    }
  }

  }
  # if no group-level paths added, don't consider z_list
  if (length(which(is.na(z_list)))==0)
  z_list <- lapply(z_list, 
                   function(x){x$sig <- ifelse(x$p < .05/n_subj, 1, 0)
                   return(x)})
  
  #subgroup based only on contemporaneous paths kmg
  if(sub_feature == "contemp"){
    mi_list <- lapply(mi_list, 
                      function(x){x[-grep('lag',mi_list[[1]]$rhs),]})
    z_list <- lapply(z_list, 
                     function(x){x[-grep('lag',z_list[[1]]$rhs),]})
  }
  
  #subgroup based only on lagged paths kmg
  if(sub_feature == "lagged"){
    mi_list <- lapply(mi_list, 
                      function(x){x[grep('lag',mi_list[[1]]$rhs),]})
    z_list <- lapply(z_list, 
                     function(x){x[grep('lag',z_list[[1]]$rhs),]})
  }
  
  # remove lines that have "NA" from z_list (occurs for exog and ar=FALSE)
  # commented out by stl april 2018 - likely to have unintended consequences
  # because it will cause off-by-one errors in the creation of the similarity matrix
  # these NA issues should instead by captured in the na.rm = T arguments
  # added to the similarity matrix. if not, we can revisit
  # z_list <- lapply(z_list, na.exclude) 
  
  sim_mi <- matrix(0, ncol = length(mi_list), nrow = length(mi_list))
  sim_z  <- sim_mi
  
  ## march 2018 stl - na.rm = TRUE added in creation of similarity matrix. 
  ## This was added to address cases where the standard errors for certain paths,
  ## particularly bidirectional paths, were NA. This caused NAs throughout
  ## the adjacency matrix. We should consider whether this is the permanent 
  ## solution we want, as it means that individuals contribute different numbers
  ## of paths to the adjacency matrix (i.e., those individuals with paths
  ## that have NA standard errors contribute fewer paths to the matrix)
  
  for (i in 1:length(mi_list)){
    for (j in 1:length(mi_list)){
      if(!hybrid && dir_prop_cutoff>0){
      sim_mi[i,j] <- sum(mi_list[[i]]$dir == 1 & mi_list[[j]]$dir == 1 & 
                           sign(mi_list[[i]]$epc) == sign(mi_list[[j]]$epc), na.rm = TRUE)
      if (length(which(is.na(z_list)))==0)
      sim_z[i,j]  <- sum(z_list[[i]]$dir == 1 & z_list[[j]]$dir == 1 &
                           sign(z_list[[i]]$z) == sign(z_list[[j]]$z), na.rm = TRUE)} else {
        sim_mi[i,j] <- sum(mi_list[[i]]$sig == 1 & mi_list[[j]]$sig == 1 & 
                             sign(mi_list[[i]]$epc) == sign(mi_list[[j]]$epc), na.rm = TRUE)
        if (length(which(is.na(z_list)))==0)
          sim_z[i,j]  <- sum(z_list[[i]]$sig == 1 & z_list[[j]]$sig == 1 &
                               sign(z_list[[i]]$z) == sign(z_list[[j]]$z), na.rm = TRUE)}
    }
  }
  
  sim           <- sim_mi + sim_z
  if (sub_sim_thresh == "search"){ # conduct search across thresholds to find lowest modularity
    res <- nloptr::nloptr(
      x0=.5,  # starting value
      m = sim,  # similarity matrix
      sub_method = sub_method,
      eval_f = modmax, # objective function
      lb = .01, # upper bound
      ub = .99, # lower bound
      opts = list(
        "algorithm"="NLOPT_LN_NELDERMEAD",
        "xtol_rel"=1.0e-8,
        maxeval = 500)
    )
    sim[which(sim <= stats::quantile(sim[upper.tri(sim, diag = FALSE)], (res$solution)))] <- 0
  } else if (sub_sim_thresh == "lowest") { # original gimme
    sim <- sim - min(sim, na.rm = TRUE)} else{
    toremove <- stats::quantile(sim[upper.tri(sim, diag = FALSE)], (sub_sim_thresh))
    sim[which(sim <= toremove)] = 0
    }
      
  diag(sim)     <- 0
  lay.sim = sim
    colnames(lay.sim) = rownames(lay.sim) = NULL
  # write.table(sim, file = paste(out_path, '/sim mat.csv', sep = ''), sep = ',')
  colnames(sim) <- rownames(sim) <- names(mi_list)
  colnames(lay.sim) <- rownames(lay.sim) <- NULL
  if(is.null(confirm_subgroup)){
    g            <- graph.adjacency(sim, mode = "undirected", weighted = TRUE)
    weights      <- E(g)$weight
    if (sub_method == "Walktrap")
      res        <- cluster_walktrap(g, weights = weights, steps = 4)
    if (sub_method == "Infomap")
      res        <- cluster_infomap(g, e.weights = weights)
    if (sub_method == "Edge Betweenness")
      res        <- cluster_edge_betweenness(g, weights = weights)
    if (sub_method == "Fast Greedy")
      res        <- cluster_fast_greedy(g, weights = weights)
    if (sub_method == "Label Prop")
      res        <- cluster_label_prop(g, weights = weights)
    if (sub_method == "Leading Eigen")
      res        <- cluster_leading_eigen(g, weights = weights)
    if (sub_method == "Louvain")
      res        <- cluster_louvain(g, weights = weights)
    if (sub_method == "Spinglass")
      res        <- cluster_spinglass(g, weights = weights)
    
    
    e = igraph::get.edgelist(igraph::graph.adjacency(lay.sim, mode = "undirected", weighted = TRUE))
    l = qgraph::qgraph.layout.fruchtermanreingold(e,
                                                  vcount = vcount(g), weights = weights/5,
                                                  area = 8*(vcount(g)^2),
                                                  repulse.rad = (vcount(g)^3.1))
    graphics::par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
    if(length(out_path)>1) {pdf(file.path(paste(out_path,'/Subgroups Plot.pdf', sep = '')))
    plot(res, g, layout = l)
    dev.off() }
    
    sub_mem    <- data.frame(names      = names(membership(res)), 
                             sub_membership = as.numeric(membership(res)))
    sub$sim         <- sim
    sub$n_subgroups <- length(unique(na.omit(sub_mem$sub_membership))) 
    sub$modularity  <- modularity(res)
    
    sub$sub_mem     <- merge(file_order, sub_mem, by = "names", all.x = TRUE)
  } else {
    sub_mem         <- confirm_subgroup
    names(sub_mem)  <- c("names", "sub_membership")
    sub$sim         <- sim
    sub$n_subgroups <- length(unique(na.omit(sub_mem$sub_membership))) 
    sub$sub_mem     <- merge(file_order, sub_mem, by = "names", all.x = TRUE)
    sub$modularity  <- modularity(graph.adjacency(sim, mode = "undirected"), (sub$sub_mem)$sub_membership)
    
    
  }
  return(sub)
}
