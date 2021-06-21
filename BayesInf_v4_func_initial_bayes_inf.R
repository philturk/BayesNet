
Inital_bayes_inf <- function(population,Network_stats, Prob_Distr, Prob_Distr_Params,
                             covPattern,
                             Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist) {
  
  G_full = graph.full(n = population, directed = FALSE)
  G_full = G_full %>% set_vertex_attr("name", value = c(0:(population-1)))
  P_start = Update_P(G_full,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist)
  
  CCMnet_Result = CCMnetpy::CCMnet_constr(Network_stats=Network_stats,
                                            Prob_Distr=Prob_Distr,
                                            Prob_Distr_Params=Prob_Distr_Params, 
                                            samplesize = as.integer(1),
                                            burnin=as.integer(100000), 
                                            interval=as.integer(10),
                                            statsonly=TRUE,
                                            G=NULL,
                                            P=NULL,
                                            population=as.integer(population), 
                                            covPattern = as.integer(covPattern),
                                            bayesian_inference = FALSE,
                                            Ia = NULL, 
                                            Il = NULL, 
                                            R = NULL, 
                                            epi_params = NULL,
                                            print_calculations = FALSE) 
  
  G_start2 = CCMnet_Result[[1]]
  
  U_P_Start = as.undirected(P_start, mode = "collapse")
  V(U_P_Start)$name <- as.character(c(0:(population-1)))
  #G_start = igraph::union(G_start2, U_P_Start)
  G_start = U_P_Start
  
  return(list(G_start, P_start))
}