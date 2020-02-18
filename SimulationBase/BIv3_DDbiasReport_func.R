DD_unbiased_sim <- function(init_seed, num_samples, biasPct) {
  outfile = paste("DD_unbiased_n", num_samples, "_", init_seed, ".RData", sep = "")
  library('igraph')   ###1/3/17
  library('CCMnet')
  library('intergraph')  ###1/3/17
  library('mvtnorm')
  library('e1071')
  
  source("BayesInf_func.R")
  
  set.seed(init_seed)
  
  Network_stats = 'DegreeDist' #'Density'
  strong_prior = TRUE
  
  beta_a = 1e-6  #run as SEIR, approx
  beta_l = 1 / 2.5
  gamma_a = 1 / 10
  gamma_l = 1 / 100
  population = 1000
  Prob_Distr_Params = list(list(rep(1, 30)))
  Prob_Distr = 'DirMult'
  bool_SIIR = TRUE
  num_init_infected = 1
  genetic_bits = 1024
  
  Prob_Distr_Params_initial = list(list(rep(1, 20)))
  
  Initial_Data = generate_epidemic_data(
    beta_a = beta_a,
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
    num_samples = num_samples
  )
  
  G = Initial_Data[[1]]
  P = Initial_Data[[2]]
  Ia = Initial_Data[[3]]
  Il = Initial_Data[[4]]
  R = Initial_Data[[5]]
  T_dist = Initial_Data[[6]]
  degtable <- Initial_Data[[7]][[1]][[1]]
  deg_samp <- rep(0:(length(degtable)-1), each = degtable)
  deg_new <- round(deg_samp*biasPct)
  samp_dist <- c(sum(deg_new == 0), tabulate(deg_new, nbins = max(max(degree(G)/2), max(deg_new))))
  Prob_Distr_Params = list(list(samp_dist+.1))
  Prob_Distr = Initial_Data[[8]]
  
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
                           covPattern = NULL,
                           remove_var_last_entry = FALSE)[[2]][[1]]
  U_graph = igraph::union(asIgraph(G_start2), as.undirected(asIgraph(P_start), mode = "collapse"))
  edgelist_Ugraph = ends(U_graph, es = c(1:ecount(U_graph)))
  G_start<-network.edgelist(edgelist_Ugraph,network.initialize(population, directed = FALSE),ignore.eval=FALSE)
  
  if (Network_stats == 'DegreeDist') {
    Prob_Distr_Params[[1]][[1]] = c(Prob_Distr_Params[[1]][[1]], rep(0.1, max(
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
  n_mcmc_trials = 1500
  n_burn = 500
  
  ecount_G = c(network.edgecount(Init_G))
  ecount_P = c(network.edgecount(Init_P))
  DD_dist = c(0)
  degdist = Init_DD
  
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
    print(mcmc_counter)
    
  }
  
  save(ecount_G, ecount_P, DD_dist, degdist, file = outfile)
}
