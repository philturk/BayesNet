
generate_epidemic_data <- function(beta_a = 1/5,
                                   beta_l = 1/2.5,
                                   gamma_a = 1/10,
                                   gamma_l = 1/100,
                                   population = 500,
                                   Prob_Distr_Params = list(list(rep(1,  10))),
                                   Prob_Distr = 'DirMult',
                                   bool_SIIR = TRUE,
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
    
    CCMnet_Result = CCMnet_constr(Network_stats=Network_stats,
                                  Prob_Distr=Prob_Distr,
                                  Prob_Distr_Params=Prob_Distr_Params, 
                                  samplesize = 2,
                                  burnin=1000000, 
                                  interval=1000,
                                  statsonly=TRUE, 
                                  P=NULL,
                                  population=population, 
                                  covPattern = covPattern,
                                  remove_var_last_entry = FALSE) 
    
    G = CCMnet_Result[[2]][[1]]
    
    Initial_Data = Initialize_G_P_Ia_Il_R(population, beta_a, beta_l, gamma_a, gamma_l, bool_SIIR, G = G, num_init_infected = num_init_infected) #This is only for testing the method
    
    G = Initial_Data[[1]]
    P = Initial_Data[[2]]
    Ia = Initial_Data[[3]]
    Il = Initial_Data[[4]]
    R = Initial_Data[[5]]
    
    if (network.edgecount(P) > 5) {
      invalid_epidemic = FALSE
    }
    
    counter = counter + 1
  }
  
  
  if (counter < max_try) {
    
    
    #####Calculate Mean and Variance of Graph properties#######
    
    if (Network_stats == 'Density') { #Assume Normal
      if (strong_prior == TRUE) {
        g_samples = sample(degree(G, gmode = 'graph'), num_samples)/(network.size(G)-1)
        Prob_Distr_Params = list(list(mean(g_samples), var(g_samples)/num_samples))
      } else {
        Prob_Distr_Params = list(list(.5, 1))
      }
      Prob_Distr = 'Normal'
    }
    
    if (Network_stats == 'Mixing') { #Assume Normal
      if (strong_prior == TRUE) {
        node_id = sample(c(1:network.size(G)), num_samples, replace = FALSE)
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
        Prob_Distr_Params[[1]][[1]] =  apply(edge_type_counts, 2, mean)
        Prob_Distr_Params[[1]][[2]] = matrix(c(apply(edge_type_counts, 2, var)[1]/num_samples,0,0,
                                               0,apply(edge_type_counts, 2, var)[2]/num_samples,0,
                                               0,0,apply(edge_type_counts, 2, var)[3]/num_samples), nrow = 3, ncol = 3) #Variance
        
      } else {
        cov_counts = tabulate(get.node.attr(G, "CovAttribute"))
        
        Prob_Distr_Params = list(NULL)
        Prob_Distr_Params[[1]][[1]] = c(.5 * (cov_counts[1]-1),
                                        .5 * (cov_counts[1]-1),
                                        .5 * (cov_counts[2]-1)) #Number or edges in mixing matrix [1,1], [1,2], and [2,2]
        Prob_Distr_Params[[1]][[2]] = matrix(c(Prob_Distr_Params[[1]][[1]][1]^2,0,0,
                                               0,Prob_Distr_Params[[1]][[1]][2]^2,0,
                                               0,0,Prob_Distr_Params[[1]][[1]][3]^2), nrow = 3, ncol = 3) #Variance
      }
      Prob_Distr = 'Normal'
    }
    
    if (Network_stats == 'DegreeDist') {  #Assume Dirichlet
      if (strong_prior == TRUE) {
        Prob_Distr_Params = list(list(tabulate(sample(degree(G, gmode = 'graph'), num_samples)+1, nbins = max(degree(G)/2)+1)))
      } else {
        Prob_Distr_Params = list(list(rep(1,  max(degree(G)/2)+1)))
      }
      Prob_Distr = 'DirMult'
    }
    
    PG_Data = Genetic_Seq_Data(P=P, Ia=Ia,final_time = max(Ia[which(Ia < Inf)]), genetic_bits = genetic_bits)
    T_dist = hamming.distance(PG_Data[,-(genetic_bits+1)])
    return(list(G, P, Ia, Il, R, T_dist, Prob_Distr_Params, Prob_Distr))
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


Initialize_G_P_Ia_Il_R <- function(population, beta_a, beta_l, gamma_a, gamma_l, bool_SIIR, G = NULL, num_init_infected = 1) {

	if (bool_SIIR == FALSE) {
		library('epinet')
		exampleepidemic <- SEIR.simulator(M = as.edgelist(G), N = population, beta = beta_l, ki = 1, ke = 1, thetae = 1/gamma_a, thetai = 1/gamma_l,latencydist='gamma')
		detach("package:epinet")		
	} else {
		exampleepidemic <- SIIR.simulator(G, population, beta_a, beta_l, gamma_a, gamma_l, num_init_infected)
	}

	P = network.initialize(population, directed = TRUE)
	P = add.edges(P, exampleepidemic[-which(is.na(exampleepidemic[,2])),2], exampleepidemic[-which(is.na(exampleepidemic[,2])),1])	

	Ia = array(Inf, dim = population)
	for (i in c(1:length(exampleepidemic[,3]))) {
		Ia[exampleepidemic[i,1]] = exampleepidemic[i,3] 
		if (bool_SIIR == FALSE) {
			Ia[exampleepidemic[i,1]] = Ia[exampleepidemic[i,1]] + abs(exampleepidemic[1,3])
			Ia[is.na(Ia)] = Inf
		}	
	}
	R = array(Inf, dim = population)
	for (i in c(1:length(exampleepidemic[,3]))) {
		R[exampleepidemic[i,1]] = exampleepidemic[i,5]
		if (bool_SIIR == FALSE) {
			R[exampleepidemic[i,1]] = R[exampleepidemic[i,1]] + abs(exampleepidemic[1,3])
			R[is.na(R)] = Inf
		}
	}	

	if (bool_SIIR == FALSE) {
		Il =  array(Inf, dim = population)
		for (i in c(1:length(Ia))) {
			if (R[i] < Inf) {
				Il[i] = runif(1,Ia[i],R[i])
				Il[is.na(Il)] = Inf
			}
		}
	} else {
		Il = array(Inf, dim = population)
		for (i in c(1:length(exampleepidemic[,3]))) {
			Il[exampleepidemic[i,1]] = exampleepidemic[i,4] #+ abs(exampleepidemic[1,3])
		}
	}

	return(list(G,P,Ia,Il,R))

}



Update_G <- function(G,P,Ia,Il,R,beta_a,beta_l, Prob_Distr_Params, Network_stats, Prob_Distr) {

  BI_Prob_Distr_Params = BI_posterior_prior(g = G,
	                                    Network_stats = Network_stats,
	                                    Prob_Distr_Params = Prob_Distr_Params) 

#   Network_stats=Network_stats
#   Prob_Distr='Normal'
#   Prob_Distr_Params=BI_Prob_Distr_Params 
#   samplesize = 2
#   burnin=100 
#   interval=2
#   statsonly=TRUE
#   TranNet = P
#   P=G
#   population=network.size(G) 
#   covPattern = NULL
#   remove_var_last_entry = FALSE
#   print_info = 0
#   BayesInference = 1
#   Ia = Ia
#   Il = Il
#   R_times = R
#   beta_a = beta_a
#   beta_l = beta_l

  covPattern = NULL
  
  if (Network_stats == 'Density') {
    BayesInference = c(1,network.edgecount(P))
  } else if (Network_stats == 'Mixing') {
    cov_counts = tabulate(get.node.attr(G, "CovAttribute"))
    BayesInference = c(1,c(summary(G~nodemix('CovAttribute'))))
    covPattern = get.node.attr(G, "CovAttribute")
  } else {
    BayesInference = c(0,0)
  }
  
  Ia[which(Ia == Inf)] = 999999+1
  Il[which(Il == Inf)] = 999999+1
  R[which(R == Inf)] = 999999+1

	CCMnet_Result = CCMnet_constr(Network_stats=Network_stats,
	                              Prob_Distr=Prob_Distr,
	                              Prob_Distr_Params=BI_Prob_Distr_Params, 
	                              samplesize = 2,
	                              burnin=500000, 
	                              interval=2,
	                              statsonly=TRUE, 
	                              P=G,
	                              population=network.size(G), 
	                              covPattern = covPattern,
	                              remove_var_last_entry = FALSE,
	                              print_info = 0,
	                              BayesInference = BayesInference,
	                              TranNet = P,
	                              Ia = Ia,
	                              Il = Il,
	                              R_times = R,
	                              beta_a = beta_a,
	                              beta_l = beta_l) 
  
	return(CCMnet_Result[[2]][[1]])
}


Update_P <- function(G,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T) {

	infected = which(R < Inf)
	infected = sample(infected, length(infected), replace = FALSE)	

	new_P = network.initialize(network.size(G), directed = TRUE)

	for (j in infected) {
		poss_parent_j = get.neighborhood(G, j, type = "combined")
		Ia_j = Ia[j]
		Il_j = Il[j]
		poss_parent_j = poss_parent_j[intersect(which(Ia[poss_parent_j] < Ia_j), which(R[poss_parent_j] > Ia_j))]
		wts = c()
		if (length(poss_parent_j) > 0) {
			for (i in poss_parent_j) {
				if (Il[i] > Ia_j) {
					w_i = beta_a*(Ia_j-Ia[i])
				} else {
					w_i = beta_l*(Ia_j-Il[i])+beta_a*(Il[i]-Ia[i])
				}
				wts = c(wts,w_i * 1/T[i,j])
			}
			parent_j_id = sample(c(1:length(poss_parent_j)),1,prob=wts)
			parent_j = poss_parent_j[parent_j_id]
	
			new_P = add.edge(new_P, tail = parent_j, head=j)

		} else {
			#orphan - or initial infected
			#print(c("ERROR: ORPHAN - ", j))
		}	
	}
	return(new_P)
}

Update_Ia <- function(G,P,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l) {

	infected = which(R < Inf)
	
	new_Ia = array(Inf,dim=length(V(G)))

	for (j in infected) {

		parent_j = unlist(neighborhood(P, 1, nodes=j-1, mode="in"))[-1] + 1
	if(length(parent_j) >0){
		parent_Ia = Ia[parent_j]
		parent_R = R[parent_j]

		contact_neighbors = unlist(neighborhood(G, 1, nodes=j-1, mode="all")) + 1
		inf_Il = Il[contact_neighbors[1]]
		neigh_Ia = Ia[contact_neighbors[-1]]
		contact_neighbors = contact_neighbors[-1][neigh_Ia< Inf]
		neigh_Ia = neigh_Ia[neigh_Ia < Inf]
		neigh_Il = Il[contact_neighbors[-1]]
		neigh_R = R[contact_neighbors[-1]]

	Ia_density <- function(x){
		acute_time = sum((neigh_Ia-x)*(neigh_Ia > x & neigh_Ia < inf_Il)) + sum((inf_Il-x)*(inf_Il < neigh_Ia)) + sum((x-neigh_Ia)*(neigh_Ia < x & x < inf_Il))
		chron_time = sum((x - neigh_Il)*(neigh_Il < x & x < neigh_R))
		
		return(exp(-gamma_a*(inf_Il-x)-beta_a*acute_time -beta_l*chron_time))
	}
		x = seq(parent_Ia, min(parent_R, inf_Il), by = 0.003)	#calculate pdf approx every day
		width = (min(parent_R, inf_Il)-parent_Ia)/length(x)
		Ia_pdf = sapply(x, Ia_density)
		Ia_pdf = Ia_pdf/(sum(Ia_pdf)*width)
		Ia_cdf = sapply(1:length(x), function(y) sum(Ia_pdf[1:y])*width)

		u = runif(1)
		new_Ia[j] = x[which(abs(Ia_cdf-u) == min(abs(Ia_cdf-u)))]
	
		}else{ new_Ia[j] = 0 }# for not just root which is by def 0...what changes when there is no parent?!
	}
	return(new_Ia)
}


Update_Il <- function(G,P,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l) {

	infected = which(R < Inf)
	
	new_Il = array(Inf,dim=length(V(G)))

	for (j in infected) {

		contact_neighbors = unlist(neighborhood(G, 1, nodes=j-1, mode="all")) + 1
		inf_Ia = Ia[contact_neighbors[1]]
		inf_R = R[contact_neighbors[1]]
		neigh_Ia = Ia[contact_neighbors[-1]]
		contact_neighbors = contact_neighbors[-1][neigh_Ia< Inf]
		neigh_Ia = neigh_Ia[neigh_Ia < Inf]
		neigh_Il = Il[contact_neighbors[-1]]
		neigh_R = R[contact_neighbors[-1]]

	Il_density <- function(x){
		acute_time = sum((x-inf_Ia)*(neigh_Ia > x & neigh_Ia < inf_R))
		chron_time = sum((neigh_Ia-x)*(neigh_Ia > x & neigh_Ia < inf_R)) + sum((inf_R-x)*(inf_R < neigh_Ia))
		
		return(exp(-gamma_a*(x-inf_Ia)-gamma_l*(inf_R-x)-beta_a*acute_time -beta_l*chron_time))
	}
		x = seq(inf_Ia, inf_R, by = 0.003)	#calculate pdf approx every day
		width = (inf_R-inf_Ia)/length(x)
		Il_pdf = sapply(x, Il_density)
		Il_pdf = Il_pdf/(sum(Il_pdf)*width)
		Il_cdf = sapply(1:length(x), function(y) sum(Il_pdf[1:y])*width)

		u = runif(1)
		new_Il[j] = x[which(abs((Il_cdf-u)) == min(abs(Il_cdf-u)))]
	
	}
	return(new_Il)
}


Update_I <- function(G,P,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l) {

	infected = which(R < Inf)
	
	new_Ia <- new_Il <- array(Inf,dim=network.size(G))

	Ia_density <- function(x){
		acute_time = sum((neigh_Ia-x)*(neigh_Ia > x & neigh_Ia < inf_Il) + (inf_Il-x)*(inf_Il < neigh_Ia) + (x-neigh_Ia)*(neigh_Ia < x & x < inf_Il))
		chron_time = sum((x - neigh_Il)*(neigh_Il < x & x < neigh_R))
		
		return(exp(-gamma_a*(inf_Il-x)-beta_a*acute_time -beta_l*chron_time))
	}

	Il_density <- function(x){
		acute_time = sum((x-inf_Ia)*(neigh_Ia > x & neigh_Ia < inf_R))
		chron_time = sum((neigh_Ia-x)*(neigh_Ia > x & neigh_Ia < inf_R)) + sum((inf_R-x)*(inf_R < neigh_Ia))
		
		return(exp(-gamma_a*(x-inf_Ia)-gamma_l*(inf_R-x)-beta_a*acute_time -beta_l*chron_time))
	}

	for (j in infected) {

		contact_neighbors = get.neighborhood(G, j, type = "combined")
		inf_Ia = Ia[j]
		inf_Il = Il[j]
		inf_R = R[j]
		neigh_Ia = Ia[contact_neighbors]
		contact_neighbors = contact_neighbors[neigh_Ia< Inf]
		neigh_Ia = neigh_Ia[neigh_Ia < Inf]
		neigh_Il = Il[contact_neighbors]
		neigh_R = R[contact_neighbors]
		parent_j = get.neighborhood(P, j, type = "in")

		if(length(parent_j) >0){
			parent_Ia = Ia[parent_j]
			parent_R = R[parent_j]
	
			x = seq(parent_Ia, min(parent_R, inf_Il), by = 0.01)
			width = (min(parent_R, inf_Il)-parent_Ia)/length(x)
			Ia_pdf = sapply(x, Ia_density)
			Ia_pdf = Ia_pdf/(sum(Ia_pdf)*width)
			Ia_cdf = sapply(1:length(x), function(y) sum(Ia_pdf[1:y])*width)
	
			u = runif(1)
			new_Ia[j] = x[which(abs(Ia_cdf-u) == min(abs(Ia_cdf-u)))[1]]
	
		}else{ 
			new_Ia[j] = 0 # for not just root which is by def 0...what changes when there is no parent?!
		}

		inf_Ia = new_Ia[j]
		x = seq(inf_Ia, inf_R, by = 0.01)	#calculate pdf approx every day
		width = (inf_R-inf_Ia)/length(x)
		Il_pdf = sapply(x, Il_density)
		Il_pdf = Il_pdf/(sum(Il_pdf)*width)
		Il_cdf = sapply(1:length(x), function(y) sum(Il_pdf[1:y])*width)

		u = runif(1)
		new_Il[j] = x[which(abs((Il_cdf-u)) == min(abs(Il_cdf-u)))[1]]
	}
	return(list(Il = new_Il, Ia = new_Ia))
}


########Phylogenetic Data Generate############33


Genetic_Seq_Data <- function(P, Ia, final_time = max(Ia[which(Ia < Inf)]), genetic_bits = 128) {

	infection_order = order(Ia)
	Genetic_Seq_mat = matrix(0, nrow = length(Ia), ncol = genetic_bits)
	update_times_a = array(0, dim = length(Ia))
	Genetic_Seq_mat = cbind(Genetic_Seq_mat, update_times_a)
	for(i in infection_order[c(1:length(which(Ia < Inf)))]) {
		parent_i = which(P[,i]==1)
		if (length(parent_i) == 0) { #root node
			Genetic_Seq_mat[i,] = c(genetic_seq_root_node(genetic_bits), Ia[i])
		} else { #has a parent
			Genetic_Seq_mat[parent_i,] = c(genetic_seq_parent_node(Genetic_Seq_mat[parent_i,-(genetic_bits+1)], Genetic_Seq_mat[parent_i,(genetic_bits+1)], Ia[i], genetic_bits), Ia[i])
			Genetic_Seq_mat[i,] = c(genetic_seq_inf_node(Genetic_Seq_mat[parent_i,-(genetic_bits+1)], Genetic_Seq_mat[parent_i,(genetic_bits+1)], Ia[i], genetic_bits), Ia[i])
		}
	}

	###Update everyone to be sampled at end of trial

	for(i in infection_order[c(1:length(which(Ia < Inf)))]) {
		Genetic_Seq_mat[i,] = c(genetic_seq_parent_node(Genetic_Seq_mat[i,-(genetic_bits+1)], Genetic_Seq_mat[i,(genetic_bits+1)], final_time, genetic_bits), final_time)
	}
	return(Genetic_Seq_mat)


}


genetic_seq_root_node <- function(genetic_bits) {
	return(sample(c(0,1), genetic_bits, replace = TRUE))
}

genetic_seq_parent_node <- function(Genetic_Seq_a, last_updatetime, ctime, genetic_bits) {
	nbit = round(ctime - last_updatetime)
	change_bits = array(0, dim = genetic_bits)
	change_bits[sample(c(1:genetic_bits), nbit, replace = FALSE)] = 1
	Genetic_Seq_a = (Genetic_Seq_a + change_bits) %% 2
	return(Genetic_Seq_a)
}

genetic_seq_inf_node <- function(Genetic_Seq_a, last_updatetime, ctime, genetic_bits) {
	nbit = 10
	change_bits = array(0, dim = genetic_bits)
	change_bits[sample(c(1:genetic_bits), nbit, replace = FALSE)] = 1
	Genetic_Seq_a = (Genetic_Seq_a + change_bits) %% 2
	return(Genetic_Seq_a)
}

#######Diagnostics###########

genetic_seq_diagnostic <- function(P,T) {

	dist_linked_cases = c()
	for (i in c(1:ecount(P))) {
		dist_linked_cases = c(dist_linked_cases, T[(P[[3]][i])+1, (P[[4]][i])+1])
	}

	d1 = density(dist_linked_cases)

	T_data = c(T[upper.tri(T,diag = FALSE)])
	d2 = density(T_data[which(T_data > 0)])

	plot(c(min(d1$x,d2$x),max(d1$x,d2$x)),c(0,max(d1$y,d2$y)), type = "n", main = "Linked vs General", xlab = "Dist in Phlyogenetics", ylab="Density")

	polygon(d1, col = '#00FF3348')
	polygon(d2, col = '#00003348')

	return(list(mean(dist_linked_cases),mean(T_data[which(T_data > 0)]), 
			var(dist_linked_cases),var(T_data[which(T_data > 0)])))
}

uncategorized_diagnostic <- function(P_a, G_a, G, P, Ia, Il, R, Initial_Data) {

	print(cor(Ia, Initial_Data[[3]]))
	print(cor(Il, Initial_Data[[4]]))

}

accuracy_G_P_diagnostic <- function(P_a, G_a, Initial_Data) {

	G = Initial_Data[[1]]
	P = Initial_Data[[2]]

	adj_m = get.adjacency(P)
	for (i in c(51:length(P_a))) {
		adj_m  = adj_m + get.adjacency(P_a[[i]])
	}
	adj_m = adj_m / length(P_a)
	adj_P = get.adjacency(P_a[[1]])

	false_positives_dist_P = adj_m[which(adj_P == 0)]
	true_positives_dist_P = adj_m[which(adj_P == 1)]

	par(mfrow = c(1,2))

	d1 = density(false_positives_dist_P, bw = .01)
	d2 = density(true_positives_dist_P, bw = .01)

	plot(c(0,1), c(0,40), type = "n", main = "Transmission Network", xlab = "Percent Edge Occurred", ylab="Density")

	polygon(d1, col = '#00FF3348')
	polygon(d2, col = '#00003348')

	print(mean(false_positives_dist_P))
	print(mean(true_positives_dist_P))

	adj_Gm = get.adjacency(G_a[[2]])
	for (i in c(3:length(G_a))) {
		adj_Gm  = adj_Gm + get.adjacency(G_a[[i]])
	}
	adj_Gm = adj_Gm / length(G_a)
	adj_G = get.adjacency(G_a[[1]])

	false_positives_dist_G = adj_Gm[which(adj_G == 0)]
	true_positives_dist_G = adj_Gm[which(adj_G == 1)]

	d1 = density(false_positives_dist_G, bw = .01)
	d2 = density(true_positives_dist_G, bw = .01)

	plot(c(0,1), c(0,40), type = "n", main = "Contact Network", xlab = "Percent Edge Occurred", ylab="Density")

	polygon(d1, col = '#00FF3348')
	polygon(d2, col = '#00003348')

	print(mean(false_positives_dist_G))
	print(mean(true_positives_dist_G))

}


network_properties_diagnostic <- function(G_a, burnin, view_graphs = FALSE) {

	burnin = burnin + 1
	population = length(network.size(G_a[[1]]))

	ecount_truth = network.edgecount(G_a[[1]])
	assort_truth =  summary(G_a[[1]]~degcor)
	degree_dist_truth = tabulate(degree(G_a[[1]])/2+1)

	ecount_a = c()
	for (i in c(burnin:length(G_a))) {
		ecount_a = c(ecount_a, network.edgecount(G_a[[i]]))
	}

	yrange = c( .99*min(c(cumsum(ecount_a) / seq_along(ecount_a), ecount_truth)) , 1.01*max(c(cumsum(ecount_a) / seq_along(ecount_a), ecount_truth)))

	if (view_graphs == TRUE) {
		plot(c(1:length(ecount_a)),cumsum(ecount_a) / seq_along(ecount_a), main = "Convergence of Density", ylab = "Cum. Mean", xlab = "Iterations", ylim = yrange )
		abline(h=ecount_truth, col = "red")
		readline(prompt = "Pause. Press <Enter> to continue...")

    plot_data = ecount_a/choose(network.size(G_a[[1]]),2)
		plot(stats::density(plot_data), main = "Density plot Network Property Edge Density", xlab = "Edge Density")
		abline(v=(ecount_truth/choose(network.size(G_a[[1]]),2)), col = "red") 
		readline(prompt = "Pause. Press <Enter> to continue...")
	}

	max_deg = 0
	for (i in c(burnin:length(G_a))) {
		if (max_deg < max(degree(G_a[[i]])/2)) {
			max_deg = max(degree(G_a[[i]])/2)
		}
	}
	max_deg = max(max_deg, max(degree(G_a[[1]])/2))

	assort_a = c()
	degree_dist_mat = matrix(0, nrow = length(c(burnin:length(G_a))), ncol = max_deg+1)
	for (i in c(burnin:length(G_a))) {
		index_i = i - burnin + 1
		assort_a = c(assort_a, summary(G_a[[i]]~degcor))
		degree_dist_mat[index_i,c(1:length(tabulate(degree(G_a[[i]])/2+1)))] = tabulate(degree(G_a[[i]])/2+1)
	}

	if (view_graphs == TRUE) {
		plot(assort_a, main = "Plot Network Property Assortativity", xlab = "Graph Index")
		abline(h=assort_truth, col = "red") 
		readline(prompt = "Pause. Press <Enter> to continue...")
	}

print(length(which(assort_a >=  assort_truth))/length(assort_a))
print(mean(assort_a))


	if (view_graphs == TRUE) {
		plot(density(assort_a), main = "Density plot Network Property Assortativity", xlab = "Assortativity Density")
		abline(v=assort_truth, col = "red") 
		readline(prompt = "Pause. Press <Enter> to continue...")
	}

####
#x = sample(assort_a,10000, replace=TRUE, prob = assort_a^10)
#mean(x)
####

#	par(mfrow = c(ceiling(sqrt(max_deg)),ceiling(sqrt(max_deg))))
#	for (i in c(1:(max_deg+1))) {
#		plot(c(1:length(degree_dist_mat[,i])),cumsum(degree_dist_mat[,i]) / seq_along(degree_dist_mat[,i]), main = "Convergence of Degree Distribution", ylab = "Cum. Mean", xlab = "Iterations")
#		abline(h=degree_dist_truth[i], col = "red") 
#	}
#	readline(prompt = "Pause. Press <Enter> to continue...")

#	par(mfrow = c(ceiling(sqrt(max_deg)),ceiling(sqrt(max_deg))))
#	for (i in c(1:(max_deg+1))) {
#		plot(density(degree_dist_mat[,i]), main = "Density plot Network Property Degree Distribution", xlab = "Degree Density")
#		abline(v=degree_dist_truth[i], col = "red") 
#	}
#	readline(prompt = "Pause. Press <Enter> to continue...")

	deg_dist_data = cbind(c(0:(max_deg)), degree_dist_truth, apply(degree_dist_mat,2,mean))
	colnames(deg_dist_data) = c("Degree", "TRUTH", "Estimated (Mean)")
#	print(deg_dist_data)
#	par(mfrow = c(1,1))
	boxplot(degree_dist_mat, main = "Degree Distribution", names = deg_dist_data[,1], xlab = "Node Degree")
	means= apply(degree_dist_mat,2,mean)
	means_true = degree_dist_truth
	points(1:length(means), means, pch = 23, cex = 0.75, bg = "blue")
	points(1:length(means_true), means_true, pch = 23, cex = 0.75, bg = "red")
	legend("topleft", legend=c("Estimated Degree Distribution (Mean)", "Truth"), pch=c(23,23), col = c("blue", "red")) 

return(list(length(which(assort_a >=  assort_truth))/length(assort_a), mean(assort_a), assort_truth, mean(ecount_a), ecount_truth))
}


BI_posterior_prior <- function(g, Network_stats, Prob_Distr_Params) {
  
  g_orig = g
  
  if  (Network_stats == 'Density') {
    
    mean_density = Prob_Distr_Params[[1]][[1]]
    var_density = Prob_Distr_Params[[1]][[2]]
    
    g_density_orig =  network.edgecount(g_orig) / choose(network.size(g_orig),2)
    var_density_orig = (g_density_orig * (1-g_density_orig))
    
    mean_density_1 = ((mean_density/var_density) + (network.edgecount(g_orig)/var_density_orig))
    mean_density_2 = ((1/var_density) + (choose(network.size(g_orig),2)/var_density_orig))
    mean_density = mean_density_1 / mean_density_2
    
    var_density = 1/mean_density_2
    
    return(list(list(mean_density, var_density)))
  }
  
  if  (Network_stats == 'Mixing') {
    
    num_cov_type = tabulate(get.node.attr(g, "CovAttribute"))
    
    mean_mixing = Prob_Distr_Params[[1]][[1]]
    var_mixing = diag(Prob_Distr_Params[[1]][[2]])
    
    mean_mixing_orig = summary(g_orig~nodemix('CovAttribute'))/c(num_cov_type[1], num_cov_type[1], num_cov_type[2])
    p1 = mean_mixing_orig[1]/(num_cov_type[1]-1)
    var_p1 = p1*(1-p1)
    
    p2 = mean_mixing_orig[2]/(num_cov_type[1]-1)
    var_p2 = p2*(1-p2)
    
    p3 = mean_mixing_orig[3]/(num_cov_type[2]-1)
    var_p3 = p3*(1-p3)
    
    var_mixing_orig = c((num_cov_type[1]-1) * var_p1,
                        (num_cov_type[1]-1) * var_p2,
                        (num_cov_type[2]-1) * var_p3)
    x_count = summary(g_orig~nodemix('CovAttribute'))
    
    mean_mixing_v = c()
    var_mixing_v = c()
    index_j = c(1,1,2)
    for (j in c(1:3)) {

      mean_density_1 = ((mean_mixing[j]/var_mixing[j]) + (x_count[j]/var_mixing_orig[j]))
      mean_density_2 = ((1/var_mixing[j]) + (num_cov_type[index_j[j]]/var_mixing_orig[j]))
      mean_mixing_v[j] = mean_density_1 / mean_density_2
    
      var_mixing_v[j] = 1/mean_density_2
    }
    
    var_mixing_mat = matrix(c(var_mixing_v[1],0,0,
                              0,var_mixing_v[2],0,
                              0,0,var_mixing_v[3]), nrow = 3, ncol = 3) #Variance
    
    return(list(list(mean_mixing_v, var_mixing_mat)))
  }
  
  if (Network_stats == 'DegreeDist') {
    
    mean_deg_dist = Prob_Distr_Params[[1]][[1]] + tabulate(degree(g_orig)/2 + 1, nbins = length(Prob_Distr_Params[[1]][[1]])) 
    return(list(list(mean_deg_dist)))
    
  }
  
  
}


#for the probability of accept of reject

update_r <- function(r, g, P, Ia, Il, R_times, edge1, beta_a, beta_l) {
  
  add_edge = 1
  if (is.adjacent(g, edge1[1], edge1[2])) { 		#Remove edge
    add_edge = 0
  } 
  
  if (add_edge == 1) {
    p_edge1 = r / (1 + r)
  } else {
    p_edge1 = 1 / (1 + r)
  }
  
  
  if ((is.adjacent(P, edge1[1], edge1[2])) || (is.adjacent(P, edge1[2], edge1[1])) ) {
    return(0)
  } else {
    if ((Ia[edge1[1]] == Inf) && (Ia[edge1[2]] == Inf)) {
      return(r)
    } else {
      Il_i = Il[edge1[1]]
      Ia_i = Ia[edge1[1]]
      R_i = R[edge1[1]]
      Il_j = Il[edge1[2]]
      Ia_j = Ia[edge1[2]]
      R_j = R[edge1[2]]
      
      if (Ia[edge1[2]] < Ia[edge1[1]]) {
        time_a = min(Il_j,Ia_i)-Ia_j
        time_l = max(min(R_j,Ia_i),Il_j) - Il_j
        muij = (1-pexp(time_a,rate = beta_a))*(1-pexp(time_l, rate = beta_l))
      } else {
        time_a = min(Il_i,Ia_j)-Ia_i
        time_l = max(min(R_i,Ia_j),Il_i) - Il_i
        muij = (1-pexp(time_a,rate = beta_a))*(1-pexp(time_l, rate = beta_l))
      }
      p_noinfect = (muij*p_edge1)/((1-p_edge1) + muij*p_edge1)
      
      if (add_edge == 1) {
        return(p_noinfect / (1 - p_noinfect) )
      } else {
        return((1-p_noinfect) / p_noinfect )
      }
    }
  }
}





