
Inital_bayes_inf <- function(population,Network_stats, Prob_Distr, Prob_Distr_Params,
                             covPattern,
                             Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist) {
  
  G_full = asNetwork(graph.full(population)) 
  P_start = Update_P(G_full,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist)
  G_start2 = CCMnet_constr(Network_stats=Network_stats,
                           Prob_Distr=Prob_Distr,
                           Prob_Distr_Params=Prob_Distr_Params,
                           samplesize = 2,
                           burnin=1000000, 
                           interval=1000,
                           statsonly=TRUE, 
                           P=NULL,
                           population=population, 
                           covPattern = covPattern,
                           remove_var_last_entry = FALSE)[[2]][[1]]
  U_graph = igraph::union(asIgraph(G_start2), as.undirected(asIgraph(P_start), mode = "collapse"))
  edgelist_Ugraph = ends(U_graph, es = c(1:ecount(U_graph)))
  G_start<-network.edgelist(edgelist_Ugraph,network.initialize(population, directed = FALSE),ignore.eval=FALSE)
  set.vertex.attribute(G_start, "CovAttribute", value = covPattern)
  
  return(list(G_start, P_start))
}