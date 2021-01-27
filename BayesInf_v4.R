
#install.packages("C:/Users/rgoyal/Desktop/Network Research/CCMnet_0.0-4.tar.gz", repos = NULL, type = "source")

library(CCMnetpy)
CCMnet_python_setup()

library('igraph')   
library('CCMnet')
library('intergraph')  

library('mvtnorm')
library('e1071')

source('BayesInf_v4_func_read_params.R')
source('BayesInf_v4_func_gen_epi_data.R')
source('BayesInf_v4_func_gen_genetic_data.R')
source('BayesInf_v4_func_initial_bayes_inf.R')
source('BayesInf_v4_func_bayes_inf.R')
source('BayesInf_v4_func_diagnostics.R')

set.seed(init_seed)

Initial_Data = generate_epidemic_data(beta_a = beta_a,
                                      beta_l = beta_l,
                                      gamma_a = gamma_a,
                                      gamma_l = gamma_l,
                                      population = population,
                                      Prob_Distr_Params = Prob_Distr_Params,
                                      Prob_Distr = Prob_Distr,
                                      num_init_infected = num_init_infected,
                                      genetic_bits = genetic_bits,
                                      Network_stats = Network_stats,
                                      strong_prior = strong_prior,
                                      num_samples = num_samples,
                                      covPattern = covPattern)

#Truth
G_truth = Initial_Data[[1]]
P_truth = Initial_Data[[2]]

#Clinical data
Ia = Initial_Data[[3]]
Il = Initial_Data[[4]]
R = Initial_Data[[5]]

#Genetic data
T_dist = Initial_Data[[6]]

#Sexual behavior data
Prob_Distr_Params = Initial_Data[[7]]
Prob_Distr = Initial_Data[[8]]

######BEGIN Bayesian Inference##################

Init_nets = Inital_bayes_inf(population,Network_stats, 
                 Prob_Distr, Prob_Distr_Params,covPattern,
                 Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist)

G = Init_nets[[1]]
P = Init_nets[[2]]

P_a = NULL
G_a = NULL

ecount_G = c()
ecount_P = c()

for (mcmc_counter in c(1:n_mcmc_trials)) {
  
  G = Update_G(G,P,Ia,Il,R,beta_a,beta_l,Prob_Distr_Params=Prob_Distr_Params, Network_stats = Network_stats, Prob_Distr = Prob_Distr)
  
  P = Update_P(G,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist)

  ecount_G = c(ecount_G, network.edgecount(G))
  ecount_P = c(ecount_P, network.edgecount(P))
  
  P_a[[mcmc_counter]] = P
  G_a[[mcmc_counter]] = G
  
  print(mcmc_counter)
  
}

plot(ecount_G)
abline(h=network.edgecount(G_truth), col = 'red')
abline(h=mean(ecount_G), col = 'blue')
