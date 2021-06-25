#SCRIPT TO SIMULATE ONE EPIDEMIC & PERFORM BAYESIAN INFERENCE
#create a Simulations directory before running (below)
#install R packages CCMnetpy, igraph, mvtnorm, e1071, tidyverse, reticulate, gtools before running (below)
#install Miniconda (run `CCMnet_python_setup()` below for prompt)
#install Python packages networkx, pandas, scipy before running (below)

####################Setup####################

#install required R packages
#remotes::install_github("ravigoyalgit/CCMnet_py")
#install.packages("igraph")
#install.packages("mvtnorm")
#install.packages("e1071")
#install.packages("tidyverse")
#install.packages("reticulate")
#install.packages("gtools")

#create directory to save simulations in
#dir.create("Simulations")

#start time
start <- Sys.time()

#pass task ID from job array to R
args <- commandArgs(trailingOnly=TRUE)
id <- as.numeric(args[[1]])

#load R packages
#library('reticulate')
library('igraph')   
library('mvtnorm')
library('e1071')
#library('tidyverse')

#install required Python packages
#py_install("networkx")
#py_install("pandas")
#py_install("scipy")

#python setup
library(CCMnetpy)
CCMnet_python_setup()

#call dependent R scripts
source('BayesInf_v4_func_read_params.R')
source('BayesInf_v4_func_gen_epi_data.R')
source('BayesInf_v4_func_gen_genetic_data.R')
source('BayesInf_v4_func_initial_bayes_inf.R')
source('BayesInf_v4_func_bayes_inf.R')
source('BayesInf_v4_func_diagnostics.R')

####################Create Network####################

set.seed(init_seed)

init_seed

save(init_seed, file=paste("Simulations/simulation-",formatC(id, width=3, format="d", flag="0"),".rda", sep=""))

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

####################Bayesian Inference####################

Init_nets = Inital_bayes_inf(population,Network_stats, 
                             Prob_Distr, Prob_Distr_Params,covPattern,
                             Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist)

G = Init_nets[[1]] #G_truth #
P = Init_nets[[2]] #P_truth #

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
  
  Prob_Distr_Params = Update_Prob_Distr_Params(G,Prob_Distr_Params_hyperprior=Prob_Distr_Params_hyperprior, Network_stats = Network_stats, Prob_Distr = Prob_Distr, Prob_Distr_Params = Prob_Distr_Params, G_stats = G_stats)
  
  print(mcmc_counter)
  
}

####################Save Results####################

#create output object
output <- list(NULL)

#list of results
results_list = list(G_stats.df,
                    Prob_Distr_Params,
                    Prob_Distr_Params_hyperprior,
                    n_mcmc_trials,
                    G_stats_truth,
                    init_seed
                    )

#add results to output
output[[3]] <- results_list

#add ID number to output
output[[1]] <- id

#end time
end <- Sys.time()

#run time
runtime <- end-start

#add runtime to output
output[[2]] <- runtime

#save simulation output with ID number
save(output, file=paste("Simulations/simulation-",formatC(id, width=3, format="d", flag="0"),".rda", sep=""))

#output (list of 3):
#[[1]] simulation ID
#[[2]] run time
#[[3]] results
