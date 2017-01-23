
deg_mixing_matrix_uni <-
  function(g, deg_max_orig) { #Statnet compliant 12/11/16
    
    deg_cor_mat <- matrix(data = 0, nrow = deg_max_orig, ncol = deg_max_orig)
    
    g_edges = unlist(g$mel)
    dim(g_edges) = c(3, length(g_edges)/3)
    
    endpoint_degree = cbind(degree(g, gmode = "graph")[g_edges[1,]], degree(g, gmode = "graph")[g_edges[2,]])
    
    for (i in c(1:network.edgecount(G))) {
      deg_cor_mat[endpoint_degree[i,1], endpoint_degree[i,2]] <- deg_cor_mat[endpoint_degree[i,1], endpoint_degree[i,2]] + 1
      if (endpoint_degree[i,2] != endpoint_degree[i,1]) {
        deg_cor_mat[endpoint_degree[i,2], endpoint_degree[i,1]] <- deg_cor_mat[endpoint_degree[i,2], endpoint_degree[i,1]] + 1
      }
    }
    
    return(deg_cor_mat)
  }
