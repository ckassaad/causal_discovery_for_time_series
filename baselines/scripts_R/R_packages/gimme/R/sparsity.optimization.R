#' @keywords internal 
modmax <- function(x, m, sub_method){
      # sparsify m based on x
      m[which(m <= stats::quantile(m[upper.tri(m, diag = FALSE)], x))] <- 0
      diag(m) <- 0
      m = graph_from_adjacency_matrix(m)
      weights      <- E(m)$weight
        if (sub_method == "Walktrap")
          modval = modularity(cluster_walktrap(m, weights = weights, steps = 4))
        if (sub_method == "Infomap")
          modval = modularity(cluster_infomap(m, e.weights = weights))
        if (sub_method == "Edge Betweenness")
          modval = modularity(cluster_edge_betweenness(m, weights = weights))
        if (sub_method == "Fast Greedy")
          modval = modularity(cluster_fast_greedy(m, weights = weights))
        if (sub_method == "Label Prop")
          modval = modularity(cluster_label_prop(m, weights = weights))
        if (sub_method == "Leading Eigen")
          modval = modularity(cluster_leading_eigen(m, weights = weights))
        if (sub_method == "Louvain")
          modval = modularity(cluster_louvain(m, weights = weights))
        if (sub_method == "Spinglass")
          modval = modularity(cluster_spinglass(m, weights = weights))
      return(-1*modval) # solvers generally want to minimze the 
}







