library(CCMnetpy)
CCMnet_python_setup()

library('igraph')   

library('mvtnorm')
library('e1071')
library('tidyverse')

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
G_stats_truth = Initial_Data[[9]]

#Clinical data
Ia = Initial_Data[[3]]
Il = Initial_Data[[4]]
R = Initial_Data[[5]]

#Genetic data
T_dist = Initial_Data[[6]]

#Hyperpriors for sexual behavior data
Prob_Distr_Params_hyperprior = Initial_Data[[7]]
Prob_Distr_Params = Initial_Data[[8]]


######BEGIN Bayesian Inference##################

Init_nets = Inital_bayes_inf(population,Network_stats, 
                             Prob_Distr, Prob_Distr_Params,covPattern,
                             Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist)

G = Init_nets[[1]] #G_truth #
P = P_truth #Init_nets[[2]] #P_truth #

P_a = NULL
G_a = NULL

G_stats.df = c()
ecount_P = c()

for (mcmc_counter in c(1:n_mcmc_trials)) {
  
  G_info = Update_G(G,P,Ia,Il,R,beta_a,beta_l,Prob_Distr_Params=Prob_Distr_Params, Network_stats = Network_stats, Prob_Distr = Prob_Distr, covPattern = covPattern)
  G = G_info[[1]]
  G_stats = as.numeric(G_info[[2]])
  
  G_stats.df = rbind(G_stats.df, G_stats)
  
  P = Update_P(G,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist)
  
  Prob_Distr_Params = Update_Prob_Distr_Params(G,Prob_Distr_Params_hyperprior=Prob_Distr_Params_hyperprior, Network_stats = Network_stats, Prob_Distr = Prob_Distr, Prob_Distr_Params = Prob_Distr_Params, G_stats = G_stats, MCMC_wgt)
  
  print(mcmc_counter)
  
}

####################

results_list = list(G_stats.df,
                    Prob_Distr_Params,
                    Prob_Distr_Params_hyperprior,
                    n_mcmc_trials,
                    G_stats_truth,
                    init_seed
                    )

