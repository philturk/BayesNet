
init_seed = runif(1,1,1000) 

#Population parameters
population = 500
covPattern = c(rep(0,population/2), rep(1,population/2)) #only use for Mixing

#Epidemic parameters
num_init_infected = 1
beta_a = 1/5
beta_l = 1/2.5
gamma_a = 1/10
gamma_l = 1/100

#Network model parameters
Network_stats = list(c("Mixing"))
Prob_Distr = list(c("Multinomial_Poisson"))
Prob_Distr_Params = vector("list", 2)
Prob_Distr_Params[[1]] = c(1000) #Number or edges in mixing matrix [1,1], [1,2], and [2,2]
Prob_Distr_Params[[2]] = c(0.3, 0.4, 0.3) 

#Prior parameters
strong_prior = TRUE
num_samples = 100

#Genetic parameters
genetic_bits = 1024

#Bayes Inf MCMC
n_mcmc_trials = 1500