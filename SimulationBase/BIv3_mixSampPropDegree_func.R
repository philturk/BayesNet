mix_prop_sim <- function(init_seed, num_samples, mixpct = 1/3) {
  outfile = paste("mix_sampProp_n", num_samples, "_m", round(mixpct*100), "_", init_seed, ".RData", sep = "")
  library('igraph')   ###1/3/17
  library('CCMnet')
  library('intergraph')  ###1/3/17
  library('mvtnorm')
  library('e1071')
  
  source("BayesInf_func.R")
  
  set.seed(init_seed)
  
  Network_stats = 'Mixing' 
  strong_prior = TRUE
  
  beta_a = 1e-6  #run as SEIR, approx
  beta_l = 1 / 2.5
  gamma_a = 1 / 10
  gamma_l = 1 / 100
  population = 1000
  
  covPattern = NULL #only use for Mixing
  
  
  Prob_Distr_Params = list(NULL)
  Prob_Distr_Params[[1]][[1]] = 20*c((1-mixpct)/2, mixpct, (1-mixpct)/2) #Number or edges in mixing matrix [1,1], [1,2], and [2,2]
  Prob_Distr_Params[[1]][[2]] = matrix(c(0.005,0,0,0,0.005,0,0,0,0.005), nrow = 3, ncol = 3) #Variance
  Prob_Distr_Params_initial = Prob_Distr_Params
  Prob_Distr = 'Normal'
  covPattern = c(rep(1,population/2), rep(2,population/2))
  
  bool_SIIR = TRUE
  num_init_infected = 1
  genetic_bits = 1024
  
  
  Initial_Data = generate_epidemic_data(beta_a = beta_a,
                                        beta_l = beta_l,
                                        gamma_a = gamma_a,
                                        gamma_l = gamma_l,
                                        population = population,
                                        Prob_Distr_Params = Prob_Distr_Params,
                                        Prob_Distr = Prob_Distr,
                                        bool_SIIR = bool_SIIR,
                                        num_init_infected = num_init_infected,
                                        genetic_bits = genetic_bits,
                                        Network_stats = Network_stats,
                                        strong_prior = strong_prior,
                                        num_samples = num_samples,
                                        covPattern = covPattern)
  
  G = Initial_Data[[1]]
  P = Initial_Data[[2]]
  Ia = Initial_Data[[3]]
  Il = Initial_Data[[4]]
  R = Initial_Data[[5]]
  T_dist = Initial_Data[[6]]
  
  degTemp <- degree(G)/ 2
  node_id = sample(c(1:network.size(G)), num_samples, replace = FALSE, prob = (degTemp+1))
  cov_types = length(unique(get.node.attr(G, "CovAttribute")))
  cov_counts = tabulate(get.node.attr(G, "CovAttribute"))
  edge_type_counts = c()
  for (i in node_id) {
    if (get.node.attr(G, "CovAttribute")[i] == 1) {
      edge_type_counts = rbind(edge_type_counts, c(tabulate(get.node.attr(G, "CovAttribute")[get.neighborhood(G, v = i, type = "combined")], nbins = cov_types),0))
    } else {
      edge_type_counts = rbind(edge_type_counts, c(0, tabulate(get.node.attr(G, "CovAttribute")[get.neighborhood(G, v = i, type = "combined")], nbins = cov_types)))
    }
  }
  svySampMean = apply(edge_type_counts, 2, mean)
  svySampVar = matrix(c(apply(edge_type_counts, 2, var)[1]/num_samples,0,0,
                        0,apply(edge_type_counts, 2, var)[2]/num_samples,0,
                        0,0,apply(edge_type_counts, 2, var)[3]/num_samples), nrow = 3, ncol = 3) #Variance
  
  Prob_Distr_Params[[1]][[1]] =  svySampMean
  Prob_Distr_Params[[1]][[2]] = svySampVar
  
  Prob_Distr = Initial_Data[[8]]
  summary(G~nodemix('CovAttribute'))
  
  ######BEGIN Bayesian Inference##################
  G_full = asNetwork(graph.full(network.size(G))) ###1/3/17
  P_start = Update_P(G_full,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist) ###1/3/17
  G_start2 = CCMnet_constr(Network_stats=Network_stats,
                           Prob_Distr=Prob_Distr,
                           Prob_Distr_Params=Prob_Distr_Params_initial,
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
  
  if (Network_stats == 'DegreeDist') {
    Prob_Distr_Params[[1]][[1]] = c(Prob_Distr_Params[[1]][[1]], rep(0.5, max(
      0, max(degree(G_start) / 2) + 1 - length(Prob_Distr_Params[[1]][[1]])
    )))
  }
  
  maxdeg = max(degree(G_start) / 2)
  
  Init_G = G
  Init_P = P
  Init_DD = summary(G~degree(0:maxdeg))
  
  G = G_start
  P = P_start
  
  P_a = list(Init_P)
  G_a = list(Init_G)
  
  mcmc_counter2 = 2 # mcmc_counter
  n_mcmc_trials = 1000
  n_burn = 50
  
  ecount_G = c(network.edgecount(Init_G))
  ecount_P = c(network.edgecount(Init_P))
  DD_dist = c(0)
  degdist = Init_DD
  mixstat = summary(Init_G~nodemix("CovAttribute"))
  
  hellinger_dist <- function(x,y){
    return(sqrt(1-sum(sqrt(x/population*y/population))))
  }
  
  for (mcmc_counter in c(mcmc_counter2:n_burn)) {
    G = Update_G(G, P, Ia, Il, R, beta_a, beta_l,
                 Prob_Distr_Params = Prob_Distr_Params,
                 Network_stats = Network_stats,
                 Prob_Distr = Prob_Distr
    )
    P = Update_P(G, Ia, Il, R, beta_a, beta_l, gamma_a, gamma_l, T_dist)
    
    print(mcmc_counter)
  }
  
  mcmc_counter2 = mcmc_counter + 1
  
  for (mcmc_counter in c(mcmc_counter2:(mcmc_counter2+n_mcmc_trials))) {
    G = Update_G(G, P, Ia, Il, R, beta_a, beta_l,
                 Prob_Distr_Params = Prob_Distr_Params,
                 Network_stats = Network_stats,
                 Prob_Distr = Prob_Distr
    )
    
    P = Update_P(G, Ia, Il, R, beta_a, beta_l, gamma_a, gamma_l, T_dist)
    
    ecount_G = c(ecount_G, network.edgecount(G))
    ecount_P = c(ecount_P, network.edgecount(P))
    DD_dist = c(DD_dist, hellinger_dist(Init_DD, summary(G~degree(0:maxdeg))))
    degdist = rbind(degdist, summary(G~degree(0:maxdeg))) 
    mixstat = rbind(mixstat, summary(G~nodemix("CovAttribute")))
    print(mcmc_counter)
    plot(ecount_G)
  }
  
  save(ecount_G, ecount_P, DD_dist, degdist, mixstat, svySampMean, svySampVar, file = outfile)
}
