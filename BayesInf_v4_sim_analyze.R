
library('tidyverse')

prior_mse_all = c()
post_mse_all = c()

for (sample_size in c(50, 100, 150, 200, 250)) {
  
  load_file = paste("C:/Users/ravij/OneDrive/Desktop/Network Research/network inference/NetBayes_git/NetBayes_git/combined-results-n", sample_size, ".rda", sep = "") 
  load(load_file)
  
  prior_mse = c()
  post_mse = c()
  
  for (i in c(1:length(results))) {
    
    G_stats.df = results[[i]][[1]]
    Prob_Distr_Params = results[[i]][[2]]
    Prob_Distr_Params_hyperprior = results[[i]][[3]]
    n_mcmc_trials = results[[i]][[4]]
    G_stats_truth = results[[i]][[5]]
    
    ProbDistr_stats.df = c()
    
    for (p_counter in c(1:n_mcmc_trials)) {
      x1 = rgamma(1, shape = Prob_Distr_Params_hyperprior[[2]][1], scale = Prob_Distr_Params_hyperprior[[2]][2])
      x2 = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
      ProbDistr_stats.df = rbind(ProbDistr_stats.df, x1*x2)
    }
    
    net_results.df = data.frame(mean_val = c(apply(G_stats.df[-c(1:800),], 2, mean), apply(ProbDistr_stats.df, 2, mean)),
                                var_val = c(apply(G_stats.df, 2, var), apply(ProbDistr_stats.df, 2, var)),
                                stat_type = c(rep("Posterior", 3), rep("Prior", 3)),
                                net_prop = rep(c("00", "01", "11"), 2),
                                truth = rep(as.numeric(G_stats_truth), 2)
    ) %>% mutate(bias2_val = (mean_val - truth)^2) %>%
      mutate(mse_val = bias2_val + var_val)
    
    # par(mfrow = c(2,2))
    # 
    # plot(G_stats.df[,1], ylab = "Mixing 0-0", xlab = "Sample",
    #      ylim = c(min(c(G_stats.df[,1], ProbDistr_stats.df[,1], as.numeric(G_stats_truth[1]))), 
    #               max(c(G_stats.df[,1], ProbDistr_stats.df[,1], as.numeric(G_stats_truth[1])))))
    # abline(h=G_stats_truth[1], col = 'red')
    # abline(h=mean(G_stats.df[,1]), col = 'blue')
    # abline(h=mean(ProbDistr_stats.df[,1]), col = 'green')
    # 
    # plot(G_stats.df[,2], ylab = "Mixing 0-1", xlab = "Sample", 
    #      ylim = c(min(c(G_stats.df[,2], ProbDistr_stats.df[,2], as.numeric(G_stats_truth[2]))), 
    #               max(c(G_stats.df[,2], ProbDistr_stats.df[,2], as.numeric(G_stats_truth[2])))))
    # abline(h=G_stats_truth[2], col = 'red')
    # abline(h=mean(G_stats.df[,2]), col = 'blue')
    # abline(h=mean(ProbDistr_stats.df[,2]), col = 'green')
    # 
    # 
    # plot(G_stats.df[,3], ylab = "Mixing 1-1", xlab = "Sample",
    #      ylim = c(min(c(G_stats.df[,3], ProbDistr_stats.df[,3], as.numeric(G_stats_truth[3]))), 
    #               max(c(G_stats.df[,3], ProbDistr_stats.df[,3], as.numeric(G_stats_truth[3])))))
    # abline(h=G_stats_truth[3], col = 'red')
    # abline(h=mean(G_stats.df[,3]), col = 'blue')
    # abline(h=mean(ProbDistr_stats.df[,3]), col = 'green')
    # 
    # par(mfrow = c(1,1))
    # 
    # 
    # ggplot(net_results.df, aes(x=net_prop, y=mse_val, fill = stat_type)) + 
    #   geom_bar(stat = "identity", position=position_dodge())
    # 
    prior_mse = c(prior_mse, net_results.df %>% filter(stat_type == "Prior") %>% pull(mse_val) %>% sum())
    post_mse = c(post_mse, net_results.df %>% filter(stat_type == "Posterior") %>% pull(mse_val) %>% sum())
    
  }
  
  prior_0 = prior_mse %>% median()
  post_0 = post_mse %>% median()
  
  #plot(prior_mse, post_mse)
  
  prior_mse_all = c(prior_mse_all, prior_0)
  post_mse_all = c(post_mse_all, post_0)
  
}

print(prior_mse_all)
print(post_mse_all)

prior_mse_all_mixing = prior_mse_all
post_mse_all_mixing = post_mse_all

################################
################################
################################
################################

prior_mse_all = c()
post_mse_all = c()

for (sample_size in c(50, 100, 150, 200, 250)) {
  
  load_file = paste("C:/Users/ravij/OneDrive/Desktop/Network Research/network inference/NetBayes_git/NetBayes_git/combined-results-n", sample_size, "-dd.rda", sep = "") 
  load(load_file)
  
  prior_mse = c()
  post_mse = c()
  
  population = 1000
  
  for (i in c(1:length(results))) {
    
    G_stats.df = results[[i]][[1]]
    Prob_Distr_Params = results[[i]][[2]]
    Prob_Distr_Params_hyperprior = results[[i]][[3]]
    n_mcmc_trials = results[[i]][[4]]
    G_stats_truth = results[[i]][[5]]
    
    ProbDistr_stats.df = c()
    
    for (p_counter in c(1:n_mcmc_trials)) {
      x2 = as.numeric(gtools::rdirichlet(n = 1, alpha = Prob_Distr_Params_hyperprior[[3]]))
      ProbDistr_stats.df = rbind(ProbDistr_stats.df, population*x2)
    }
    
    
    net_results_full.df = data.frame(mean_val = c(apply(G_stats.df[c(801:1000),], 2, mean), apply(ProbDistr_stats.df, 2, mean)),
                                     var_val = c(apply(G_stats.df, 2, var), apply(ProbDistr_stats.df, 2, var)),
                                     stat_type = c(rep("Posterior", population), rep("Prior", population)),
                                     net_prop = rep(paste("Degree", c(0:(population-1)), sep="_"), 2),
                                     truth = rep(as.numeric(G_stats_truth), 2)
    ) %>% mutate(bias2_val = (mean_val - truth)^2) %>%
      mutate(mse_val = bias2_val + var_val)
    
    net_results.df = net_results_full.df %>% filter(net_prop %in% paste("Degree", c(2:7), sep="_"))
    
    prior_mse = c(prior_mse, net_results.df %>% filter(stat_type == "Prior") %>% pull(mse_val) %>% sum())
    post_mse = c(post_mse, net_results.df %>% filter(stat_type == "Posterior") %>% pull(mse_val) %>% sum())
    
    print(i)
  }
  
  post_mse %>% hist()
  
  prior_0 = prior_mse %>% mean()
  post_0 = post_mse %>% mean()
  
  prior_mse_all = c(prior_mse_all, prior_0)
  post_mse_all = c(post_mse_all, post_0)
  
}

print(prior_mse_all)
print(post_mse_all)

prior_mse_all_degree = prior_mse_all
post_mse_all_degree = post_mse_all