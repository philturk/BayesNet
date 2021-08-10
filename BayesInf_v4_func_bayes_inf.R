Update_G <- function(G,P,Ia,Il,R,beta_a,beta_l, Prob_Distr_Params, Network_stats, Prob_Distr, covPattern) {
  
  Ia[which(Ia == Inf)] = 999999+1
  Il[which(Il == Inf)] = 999999+1
  R[which(R == Inf)] = 999999+1
  
  G_df = igraph::as_data_frame(G)
  G_df[,1] =  as.integer(G_df[,1])
  G_df[,2] =  as.integer(G_df[,2])
  
  P_df = igraph::as_data_frame(P)
  P_df[,1] =  as.integer(P_df[,1])
  P_df[,2] =  as.integer(P_df[,2])
  
  epi_params = c(beta_a, beta_l)
  
  CCMnet_Result = CCMnetpy::CCMnet_constr(Network_stats=Network_stats,
                                          Prob_Distr=Prob_Distr,
                                          Prob_Distr_Params=Prob_Distr_Params, 
                                          samplesize = as.integer(1),
                                          burnin=as.integer(10000), 
                                          interval=as.integer(10),
                                          statsonly=TRUE,
                                          G=G_df,
                                          P=P_df,
                                          population=as.integer(population), 
                                          covPattern = as.integer(covPattern),
                                          bayesian_inference = TRUE,
                                          Ia = Ia, 
                                          Il = Il, 
                                          R = R, 
                                          epi_params = epi_params,
                                          print_calculations = FALSE) 
  
  return(list(CCMnet_Result[[1]], CCMnet_Result[[2]]))
}


Update_P <- function(G,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T) {
  
  infected = which(R < Inf)
  infected = sample(infected, length(infected), replace = FALSE)	
  
  new_P = data.frame(from = NULL, to = NULL)
  
  for (j in infected) {
    poss_parent_j = igraph::neighbors(G, j) #takes vertex id, returns names = vertex_id - 1
    poss_parent_j = as.numeric(V(G)[poss_parent_j])
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
      
      new_P_edge = data.frame(from = parent_j-1, to = j-1)
      new_P = rbind(new_P, new_P_edge)
    } else {
      #orphan - or initial infected
      #print(c("ERROR: ORPHAN - ", j))
    }	
  }
  
  nodes_attr_df = data.frame(name = c(0:(population-1)))
  P = graph_from_data_frame(new_P, directed=TRUE, vertices = nodes_attr_df)
  
  return(P)
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



Update_Prob_Distr_Params <- function(g, Prob_Distr_Params_hyperprior, Network_stats, Prob_Distr, Prob_Distr_Params, G_stats, MCMC_wgt) {
  
  if  (Network_stats == 'Mixing') {
    if (Prob_Distr[[1]][1] == 'Multinomial_Poisson') {
      if (Prob_Distr_Params_hyperprior[[1]][1] == 'Dirichlet_Gamma') {
        gamma_kappa = Prob_Distr_Params_hyperprior[[2]][1]
        gamma_theta = Prob_Distr_Params_hyperprior[[2]][2]
        
        alpha = Prob_Distr_Params_hyperprior[[3]]
        
        g_edgecount = igraph::gsize(g)
        
        gamma_kappa_update = gamma_kappa + g_edgecount
        gamma_theta_update =  gamma_theta/(1*gamma_theta + 1)
        
        alpha_update = alpha + (G_stats * MCMC_wgt)

        Prob_Distr_Params[[1]][1] = rgamma(1, shape = gamma_kappa_update, scale = gamma_theta_update)
        Prob_Distr_Params[[2]] = as.numeric(gtools::rdirichlet(n = 1, alpha = alpha_update))
      }
    }
  }
  
  if  (Network_stats == 'Degree') {
    if (Prob_Distr[[1]][1] == 'Multinomial_Poisson') {
      if (Prob_Distr_Params_hyperprior[[1]][1] == 'Dirichlet_Gamma') {

        alpha = Prob_Distr_Params_hyperprior[[3]]
        
        v_deg =igraph::degree(g, mode = "all")
        g_deghist = ((v_deg) %>% tabulate(nbins = population)) * MCMC_wgt
        
        Prob_Distr_Params[[2]] = as.numeric(gtools::rdirichlet(n = 1, alpha = alpha + g_deghist))
        min_val = min(Prob_Distr_Params[[2]][which(Prob_Distr_Params[[2]]>0)])
        Prob_Distr_Params[[2]][which(Prob_Distr_Params[[2]]==0)] = min_val
        Prob_Distr_Params[[2]] = Prob_Distr_Params[[2]]/sum(Prob_Distr_Params[[2]])

      }
    }
  }
  
  return(Prob_Distr_Params)
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

