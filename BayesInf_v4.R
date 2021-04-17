
#install.packages("C:/Users/rgoyal/Desktop/Network Research/CCMnet_0.0-4.tar.gz", repos = NULL, type = "source")

library(CCMnetpy)
CCMnet_python_setup()

library('igraph')   
#library('CCMnet')
#library('intergraph')  

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
  
  #ecount_G = c(ecount_G, igraph::ecount(G))
  ecount_P = c(ecount_P, igraph::ecount(P))
  
  #P_a[[mcmc_counter]] = P
  #G_a[[mcmc_counter]] = G
  
  print(mcmc_counter)
  
}

####################

ProbDistr_stats.df = c()

for (p_counter in c(1:n_mcmc_trials)) {
  x1 = rgamma(1, shape = Prob_Distr_Params_hyperprior[[2]][1], scale = Prob_Distr_Params_hyperprior[[2]][2])
  x2 = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
  ProbDistr_stats.df = rbind(ProbDistr_stats.df, x1*x2)
  print(p_counter)
}

par(mfrow = c(2,2))

plot(G_stats.df[,1], ylab = "Mixing 0-0", xlab = "Sample",
     ylim = c(min(c(G_stats.df[,1], ProbDistr_stats.df[,1], as.numeric(G_stats_truth[1]))), 
              max(c(G_stats.df[,1], ProbDistr_stats.df[,1], as.numeric(G_stats_truth[1])))))
abline(h=G_stats_truth[1], col = 'red')
abline(h=mean(G_stats.df[,1]), col = 'blue')
#abline(h=mean(G_stats.df[,1]) + sd(G_stats.df[,1]), col = 'blue', lty = 2)
#abline(h=mean(G_stats.df[,1]) - sd(G_stats.df[,1]), col = 'blue', lty = 2)
abline(h=mean(ProbDistr_stats.df[,1]), col = 'green')
#abline(h=mean(ProbDistr_stats.df[,1]) + sd(ProbDistr_stats.df[,1]), col = 'green', lty = 2)
#abline(h=mean(ProbDistr_stats.df[,1]) - sd(ProbDistr_stats.df[,1]), col = 'green', lty = 2)


plot(G_stats.df[,2], ylab = "Mixing 0-1", xlab = "Sample", 
     ylim = c(min(c(G_stats.df[,2], ProbDistr_stats.df[,2], as.numeric(G_stats_truth[2]))), 
                                                                    max(c(G_stats.df[,2], ProbDistr_stats.df[,2], as.numeric(G_stats_truth[2])))))
abline(h=G_stats_truth[2], col = 'red')
abline(h=mean(G_stats.df[,2]), col = 'blue')
#abline(h=mean(G_stats.df[,2]) + sd(G_stats.df[,2]), col = 'blue', lty = 2)
#abline(h=mean(G_stats.df[,2]) - sd(G_stats.df[,2]), col = 'blue', lty = 2)
abline(h=mean(ProbDistr_stats.df[,2]), col = 'green')
#abline(h=mean(ProbDistr_stats.df[,2]) + sd(ProbDistr_stats.df[,2]), col = 'green', lty = 2)
#abline(h=mean(ProbDistr_stats.df[,2]) - sd(ProbDistr_stats.df[,2]), col = 'green', lty = 2)

plot(G_stats.df[,3], ylab = "Mixing 1-1", xlab = "Sample",
     ylim = c(min(c(G_stats.df[,3], ProbDistr_stats.df[,3], as.numeric(G_stats_truth[3]))), 
              max(c(G_stats.df[,3], ProbDistr_stats.df[,3], as.numeric(G_stats_truth[3])))))
abline(h=G_stats_truth[3], col = 'red')
abline(h=mean(G_stats.df[,3]), col = 'blue')
#abline(h=mean(G_stats.df[,3]) + sd(G_stats.df[,3]), col = 'blue', lty = 2)
#abline(h=mean(G_stats.df[,3]) - sd(G_stats.df[,3]), col = 'blue', lty = 2)
abline(h=mean(ProbDistr_stats.df[,3]), col = 'green')
#abline(h=mean(ProbDistr_stats.df[,3]) + sd(ProbDistr_stats.df[,3]), col = 'green', lty = 2)
#abline(h=mean(ProbDistr_stats.df[,3]) - sd(ProbDistr_stats.df[,3]), col = 'green', lty = 2)

####################

library('ggplot2')

net_results.df = data.frame(mean_val = c(apply(G_stats.df, 2, mean), apply(ProbDistr_stats.df, 2, mean)),
           var_val = c(apply(G_stats.df, 2, var), apply(ProbDistr_stats.df, 2, var)),
           stat_type = c(rep("Posterior", 3), rep("Prior", 3)),
           net_prop = rep(c("00", "01", "11"), 2),
           truth = rep(as.numeric(G_stats_truth), 2)
) %>% mutate(bias2_val = (mean_val - truth)^2) %>%
  mutate(mse_val = bias2_val + var_val)

ggplot(net_results.df, aes(x=net_prop, y=mean_val, fill = stat_type)) + 
  geom_bar(stat = "identity", position=position_dodge())

ggplot(net_results.df, aes(x=net_prop, y=var_val, fill = stat_type)) + 
  geom_bar(stat = "identity", position=position_dodge())

ggplot(net_results.df, aes(x=net_prop, y=bias2_val, fill = stat_type)) + 
  geom_bar(stat = "identity", position=position_dodge())

ggplot(net_results.df, aes(x=net_prop, y=mse_val, fill = stat_type)) + 
  geom_bar(stat = "identity", position=position_dodge())


#################################

cbind(((apply(G_stats.df, 2, mean) %>% c())/population),
      (as.numeric(G_stats_truth[1,]))/population) %>% View()

ProbDistr_stats.df = c()
for (p_counter in c(1:n_mcmc_trials)) {
  x2 = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
  ProbDistr_stats.df = rbind(ProbDistr_stats.df, population*x2)
}

par(mfrow = c(3,3))

for (i in c(1:9)) {
  plot(G_stats.df[,i], ylab = paste("Degree", i), xlab = "Sample",
       ylim = c(min(c(G_stats.df[,i], ProbDistr_stats.df[,i], as.numeric(G_stats_truth[i]))), 
                max(c(G_stats.df[,i], ProbDistr_stats.df[,i], as.numeric(G_stats_truth[i])))))
  abline(h=G_stats_truth[i], col = 'red')
  abline(h=mean(G_stats.df[,i]), col = 'blue')
  abline(h=mean(ProbDistr_stats.df[,i]), col = 'green')
}

net_results_full.df = data.frame(mean_val = c(apply(G_stats.df, 2, mean), apply(ProbDistr_stats.df, 2, mean)),
                            var_val = c(apply(G_stats.df, 2, var), apply(ProbDistr_stats.df, 2, var)),
                            stat_type = c(rep("Posterior", population), rep("Prior", population)),
                            net_prop = rep(paste("Degree", c(0:99), sep="_"), 2),
                            truth = rep(as.numeric(G_stats_truth), 2)
) %>% mutate(bias2_val = (mean_val - truth)^2) %>%
  mutate(mse_val = bias2_val + var_val)

net_results.df = net_results_full.df %>% filter(net_prop %in% paste("Degree", c(2:7), sep="_"))

ggplot(net_results.df, aes(x=net_prop, y=mean_val, fill = stat_type)) + 
  geom_bar(stat = "identity", position=position_dodge())

ggplot(net_results.df, aes(x=net_prop, y=var_val, fill = stat_type)) + 
  geom_bar(stat = "identity", position=position_dodge())

ggplot(net_results.df, aes(x=net_prop, y=bias2_val, fill = stat_type)) + 
  geom_bar(stat = "identity", position=position_dodge())

ggplot(net_results.df, aes(x=net_prop, y=mse_val, fill = stat_type)) + 
  geom_bar(stat = "identity", position=position_dodge())


##################################
