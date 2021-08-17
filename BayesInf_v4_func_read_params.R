
init_seed = runif(1,1,1000) 

#Population parameters
population = 1000
covPattern = c(rep(0,population/2), rep(1,population/2)) #only use for Mixing

#Epidemic parameters
num_init_infected = 1
beta_a = 0
beta_l = 1
gamma_a = 1/5
gamma_l = 1/2

#Network model parameters
#Network_stats = list(c("Mixing"))
Network_stats = list(c("Degree"))
Prob_Distr = list(c("Multinomial_Poisson"))
Prob_Distr_Params = vector("list", 2)
Prob_Distr_Params[[1]] = c(1000) #Number or edges in mixing matrix [1,1], [1,2], and [2,2]
#Prob_Distr_Params[[2]] = c(0.1, 0.3, 0.6) 
Prob_Distr_Params[[2]] = c(c(0.0001, 0.0001, 0.2, 0.3, 0.2, 0.1, 0.05, 0.025, 0.025), rep(0.0001, 491)) 


#Hyper prior parameters
Prob_Distr_Params_hyperprior = vector("list", 3)
Prob_Distr_Params_hyperprior[[1]][1] = 'Dirichlet_Gamma'
Prob_Distr_Params_hyperprior[[2]][1] = 1 #gamma k (shape)
Prob_Distr_Params_hyperprior[[2]][2] = 1000 #gamma theta (scale)
#Prob_Distr_Params_hyperprior[[3]] = rep(1/length(Prob_Distr_Params[[2]]), length(Prob_Distr_Params[[2]])) #dirichlet alpha
Prob_Distr_Params_hyperprior[[3]] = rep(1/population, population) #dirichlet alpha


strong_prior = TRUE
num_samples = 100

#Genetic parameters
genetic_bits = 1024

#Bayes Inf MCMC
n_mcmc_trials = 1000

#Prior vs. MCMC weight
MCMC_wgt = 0

#set transmission network to truth
init_P_truth_bool = TRUE
