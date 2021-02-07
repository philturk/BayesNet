generate_epidemic_data <- function(beta_a = 1/5,
                                   beta_l = 1/2.5,
                                   gamma_a = 1/10,
                                   gamma_l = 1/100,
                                   population = 500,
                                   Prob_Distr_Params = list(list(rep(1,  10))),
                                   Prob_Distr = 'DirMult',
                                   num_init_infected = 1,
                                   genetic_bits = 1024,
                                   Network_stats = NULL,
                                   strong_prior = NULL,
                                   num_samples = 100,
                                   covPattern = NULL) {
  
  invalid_epidemic = TRUE
  counter = 1
  max_try = 5
  
  while(invalid_epidemic & counter < max_try) {

    CCMnet_Result = CCMnetpy::CCMnet_constr(Network_stats=Network_stats,
                                  Prob_Distr=Prob_Distr,
                                  Prob_Distr_Params=Prob_Distr_Params, 
                                  samplesize = as.integer(1),
                                  burnin=as.integer(1000000), 
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
    
    G = CCMnet_Result[[1]]
    
    Initial_Data = Initialize_G_P_Ia_Il_R(population, beta_a, beta_l, gamma_a, gamma_l, G = G, num_init_infected = num_init_infected) 
    
    G = Initial_Data[[1]]
    P = Initial_Data[[2]]
    Ia = Initial_Data[[3]]
    Il = Initial_Data[[4]]
    R = Initial_Data[[5]]
    
    if (igraph::ecount(P) > 5) {
      invalid_epidemic = FALSE
    }
    
    counter = counter + 1
  }
  
  
  if (counter < max_try) {
  
    #####Calculate Mean and Variance of Graph properties#######
    
    if (Network_stats == 'Mixing') { 
      if (strong_prior == TRUE) {
        gamma_kappa = Prob_Distr_Params_hyperprior[[2]][1]
        gamma_theta = Prob_Distr_Params_hyperprior[[2]][2]
        
        alpha = Prob_Distr_Params_hyperprior[[3]]
        
        g_edgecount = igraph::gsize(G)
        
        Prob_Distr_Params_hyperprior[[2]][1] = gamma_kappa + g_edgecount
        Prob_Distr_Params_hyperprior[[2]][2] = gamma_theta/(1*gamma_theta + 1)
        
        Prob_Distr_Params_hyperprior[[3]] = alpha + as.numeric(CCMnet_Result[[2]])
        
      }
      Prob_Distr_Params[[1]][1] = rgamma(1, shape = Prob_Distr_Params_hyperprior[[2]][1], scale = Prob_Distr_Params_hyperprior[[2]][2])
      Prob_Distr_Params[[2]] = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
    }
    
    PG_Data = Genetic_Seq_Data(P=P, Ia=Ia,final_time = max(Ia[which(Ia < Inf)]), genetic_bits = genetic_bits)
    T_dist = hamming.distance(PG_Data[,-(genetic_bits+1)])
    
    return(list(G, P, Ia, Il, R, T_dist, Prob_Distr_Params_hyperprior, Prob_Distr_Params))
  } else {
    print("ERROR: FAILED")
    return(list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL))
  }
}

SIIR.simulator <- function(g, population, beta_a, beta_l, gamma_a, gamma_l, num_init_infected = 1) {
  init_infected = sample(c(1:population),num_init_infected,replace=FALSE)
  
  Ia = array(Inf, dim=population)
  Il = array(Inf, dim=population)
  R = array(Inf, dim=population)
  Parent = array(NA, dim=population)
  
  for (init in init_infected) {
    Ia[init] = 0
    Il[init] = rexp(1, gamma_a)
    R[init] = Il[init] + rexp(1, gamma_l)
    Parent[init] = NA
  }
  
  Infected_Nodes = init_infected
  while(length(Infected_Nodes) > 0) {
    
    Node_index = min(which(Ia[Infected_Nodes] == min(Ia[Infected_Nodes]))) #identify and remove more current infected
    Node_id = Infected_Nodes[Node_index]
    Infected_Nodes = Infected_Nodes[-Node_index]
    
    Node_neighbors = which(g[Node_id,]==1)
    if (length(Node_neighbors) > 0) {
      for (poss_infect in Node_neighbors) {
        acute_time = rexp(1, beta_a) #time of infectious contact during acute
        if ( ((acute_time + Ia[Node_id]) < Il[Node_id])  && ((acute_time + Ia[Node_id]) < Ia[poss_infect]) )  { #still in acute phase
          Ia[poss_infect] = acute_time + Ia[Node_id]
          Il[poss_infect] = rexp(1, gamma_a) + Ia[poss_infect]
          R[poss_infect] = Il[poss_infect] + rexp(1, gamma_l)
          Parent[poss_infect] = Node_id
          Infected_Nodes = unique(c(Infected_Nodes, poss_infect))
        } else {
          longterm_time = rexp(1, beta_l) #time of infectious contact during long-term
          if ( ((longterm_time + Il[Node_id]) < R[Node_id]) && ((longterm_time + Il[Node_id]) < Ia[poss_infect]) ) { #still in long-term phase
            Ia[poss_infect] = longterm_time + Il[Node_id]
            Il[poss_infect] = rexp(1, gamma_a) + Ia[poss_infect]
            R[poss_infect] = Il[poss_infect] + rexp(1, gamma_l)
            Parent[poss_infect] = Node_id
            Infected_Nodes = unique(c(Infected_Nodes, poss_infect))
          }
        }
      }
    }
  }
  infected_nodes = which(Ia < Inf)
  exampleepidemic = matrix(nrow = length(infected_nodes), ncol = 5)
  counter = 1
  for (infected_node in infected_nodes) {
    exampleepidemic[counter, ] = c(infected_node, Parent[infected_node], Ia[infected_node], Il[infected_node], R[infected_node])
    counter = counter + 1
  }
  colnames(exampleepidemic) = c("Infected", "Parent", "Acute", "Long-Term","Trt")
  return(exampleepidemic)
}


Initialize_G_P_Ia_Il_R <- function(population, beta_a, beta_l, gamma_a, gamma_l, G = NULL, num_init_infected = 1) {
  
  exampleepidemic <- SIIR.simulator(G, population, beta_a, beta_l, gamma_a, gamma_l, num_init_infected)
  
  P = igraph::make_empty_graph(n = population, directed = TRUE)
  
  P_edges = c(rbind(exampleepidemic[-which(is.na(exampleepidemic[,2])),2], exampleepidemic[-which(is.na(exampleepidemic[,2])),1]))
  P = igraph::add_edges(P, P_edges)
  
  Ia = array(Inf, dim = population)
  for (i in c(1:length(exampleepidemic[,3]))) {
    Ia[exampleepidemic[i,1]] = exampleepidemic[i,3] 
  }
  
  R = array(Inf, dim = population)
  for (i in c(1:length(exampleepidemic[,3]))) {
    R[exampleepidemic[i,1]] = exampleepidemic[i,5]
  }	
  
  Il = array(Inf, dim = population)
  for (i in c(1:length(exampleepidemic[,3]))) {
    Il[exampleepidemic[i,1]] = exampleepidemic[i,4] #+ abs(exampleepidemic[1,3])
  }

  
  return(list(G,P,Ia,Il,R))
  
}