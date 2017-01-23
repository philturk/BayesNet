
#install.packages("C:/Users/rgoyal/Desktop/Network Research/CCMnet_0.0-4.tar.gz", repos = NULL, type = "source")

library('igraph')   ###1/3/17
library('CCMnet')
library('intergraph')  ###1/3/17

setwd('C:\\Users\\rgoyal\\Desktop\\Network Research\\network inference\\BayesNet-master')

source("Generate_Network_Data.R")
source("BayesInf_func_C_3.R")
source("Sample_Network_Data_func.R")
source("GUF_BayesInf_func.R")

source("deg_mixing_matrix_uni.R")  ###1/3/17


library('mvtnorm')
library('e1071')

#Fixed Values

beta_a = 1/5
beta_l = 1/2.5
gamma_a = 1/10
gamma_l = 1/100

population = 500

###########Initial G###########

bool_ER = TRUE
bool_Assort = FALSE
bool_Clustering = FALSE

ER_prob = .005
Assort_val = .2
Cluster_val = 5

mean_bias = c()

#for (counter_ex in c(1:10)) {

G = G.generate(g_type = 2, population = population, bool_ER = bool_ER, bool_Assort = bool_Assort, 
               bool_Clustering = bool_Clustering, ER_prob = ER_prob, Assort_val = Assort_val, Cluster_val = Cluster_val) 

summary(G~triangle)
tabulate(degree(G)/2 + 1)
network.edgecount(G)
summary(G~degcor)

#######Initial Data: (P,Ia,Il, R)#########

bool_SIIR = TRUE
num_init_infected = 1

init_seed = runif(1,1,1000) #1.100466 #757.2034 #546.148 #633.2773 # 563.1602 #

Initial_Data = Initialize_G_P_Ia_Il_R(init_seed, population, beta_a, beta_l, gamma_a, gamma_l, bool_SIIR, ER_prob, bool_Clustering = bool_Clustering, Cluster_val = Cluster_val, G = G, num_init_infected = num_init_infected) #This is only for testing the method

G = Initial_Data[[1]]
P = Initial_Data[[2]]
Ia = Initial_Data[[3]]
Il = Initial_Data[[4]]
R = Initial_Data[[5]]

network.edgecount(P)


#######Calcuate Genetic Sequence Data########

genetic_bits = 1024

if (network.edgecount(P) > 5) {
  PG_Data = Genetic_Seq_Data(P=P, Ia=Ia,final_time = max(Ia[which(Ia < Inf)]), genetic_bits = genetic_bits)
  T_dist = hamming.distance(PG_Data[,-(genetic_bits+1)])
}

#####Calculate Mean and Variance of Graph properties#######

sample_deg_dist = FALSE #TRUE #
sample_density = FALSE #FALSE #
sample_dmm = TRUE #FALSE #
num_samples = 100

strong_prior = TRUE

eta = list(NULL,NULL,NULL,NULL,NULL,NULL)

##Should debug these call to be sure getting right values - GUF fns use igraph commands
if (sample_density) {
  
  
  eta[[1]] = .5 #network parameter and uncertainty
  eta[[2]] = 1 
  
  if (strong_prior == TRUE) {
    eta[[1]] = network.edgecount(G)/choose(population,2) #network parameter and uncertainty
    eta[[2]] = .0001 
  }
  Prob_Distr_Params = list(list(eta[[1]], eta[[2]]))
  Network_stats = 'Density'
}
if (sample_deg_dist) {

  eta[[3]] = population * array(rep(1/(max(degree(G)/2)+1),  max(degree(G)/2)+1 ), dim = max(degree(G)/2)+1)
  
  eta[[4]] = population * var_deg_dist_compute(eta[[3]])
  
  if (strong_prior == TRUE) {
    eta[[3]] = tabulate(degree(G)/2+1) + .25 #MODIFIED
    eta[[4]] = (population * var_deg_dist_compute(eta[[3]])) / num_samples
  }
  Prob_Distr_Params = list(list(eta[[3]], eta[[4]]))
  Network_stats = 'DegreeDist'
}
if (sample_dmm) {

  eta[[5]] = matrix(data = rep(1/(max(degree(G)/2) * (max(degree(G)/2)+1)), (max(degree(G)/2) * (max(degree(G)/2)))), nrow = max(degree(G)/2), ncol = (max(degree(G)/2)))
  eta[[6]] = var_deg_dist_compute(eta[[5]][upper.tri(eta[[5]], diag = TRUE)])
  
  if (strong_prior == TRUE) {
    eta[[5]] = deg_mixing_matrix_uni(G, max(degree(G)/2)) ###1/3/17
    
    rand_eta = runif(n = length(c(eta[[5]])), min = 0, max = 1)
    eta[[5]] = eta[[5]] + rand_eta  ###1/22/17
    eta[[5]] = eta[[5]] / sum(eta[[5]][upper.tri(eta[[5]], diag = TRUE)]) #network.edgecount(G)  ###1/22/17
    
    eta[[6]] = var_deg_dist_compute(eta[[5]][upper.tri(eta[[5]], diag = TRUE)]) / num_samples ###1/22/17
    
  }
  Prob_Distr_Params = list(list(eta[[5]], eta[[6]]))
  Network_stats = 'DegMixing'
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#statnet compatible to here
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


######BEGIN Bayesian Inference##################

valid_graph_start = FALSE
G_full = asNetwork(graph.full(network.size(G))) ###1/3/17
counter = 1  ###1/3/17
while(valid_graph_start == FALSE) {
  G_start = G.generate(g_type = 2, population = population, bool_ER = TRUE, bool_Assort = FALSE, 
                       bool_Clustering = FALSE, ER_prob = ER_prob, Assort_val = Assort_val, Cluster_val = Cluster_val) ###1/3/17
  P_start = Update_P(G_full,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist) ###1/3/17
  G_start = asNetwork(union(asIgraph(G_start), as.undirected(asIgraph(P_start), mode = "collapse"))) ###1/3/17
  if (max(degree(G_start)/2) <= max(degree(G)/2)) {
    valid_graph_start = TRUE
  } else { ###1/3/17
    print(paste("Failed Trial: ", counter)) ###1/3/17
    counter = counter + 1 ###1/3/17
    ER_prob = ER_prob / 2  ###1/3/17
  } ###1/3/17
}


Init_G = G
Init_P = P

G = G_start
P = P_start

P_a = list(Init_P)
G_a = list(Init_G)

mcmc_counter2 = 2 # mcmc_counter 
n_mcmc_trials = 1500

ecount_G = c(network.edgecount(Init_G))
ecount_P = c(network.edgecount(Init_P))

for (mcmc_counter in c(mcmc_counter2:n_mcmc_trials)) {
  
  G = Update_G(G,P,Ia,Il,R,beta_a,beta_l,Prob_Distr_Params=Prob_Distr_Params, Network_stats = Network_stats)
  
  P = Update_P(G,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l, T_dist)
  
  #  Ia = Update_Ia(G,P,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l)	
  #	Il = Update_Il(G,P,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l)	
  #	new_I = Update_I(G,P,Ia,Il,R,beta_a,beta_l,gamma_a,gamma_l)
  #	Ia = new_I$Ia
  #	Il = new_I$Il
  
  ecount_G = c(ecount_G, network.edgecount(G))
  ecount_P = c(ecount_P, network.edgecount(P))
  
  P_a[[mcmc_counter]] = P
  G_a[[mcmc_counter]] = G
  
  print(mcmc_counter)
  
}

mean_bias = c(mean_bias, ecount_G[1] - mean(ecount_G[c(1000:1500)]))



plot(ecount_G)
abline(h=ecount_G[1], col = 'red')
abline(h=mean(ecount_G[-1]), col = 'blue')

plot(ecount_P)
abline(h=ecount_P[1], col = 'red')
abline(h=mean(ecount_P[-1]), col = 'blue')


prob_type_s = paste(prob_type[1],prob_type[2],prob_type[3],prob_type[4],prob_type[5], sep = "")
file_timestamp = round(as.numeric(Sys.time()))

if (run_program_option <= 1) {
  result_file = paste("BayesNetInf_",population,"_",ER_prob,"_",Assort_val,"_",prob_type_s,"_",file_timestamp,".RData", sep="")
}
if (run_program_option == 2) {
  result_file = paste("/home/rg88/BayesianNetworkInference/results/BayesNetInf_",population,"_",ER_prob,"_",Assort_val,"_",prob_type_s,"_",file_timestamp,".RData", sep="")
}

save(P_a, G_a, Initial_Data, file = result_file)

file_dir = "1_21_13"
file_dir = "2_06_13"
file_dir = "Apr122013"
file_dir = "Apr142013"
result_file_a <- list.files(paste("C:/Users/Ravi/Desktop/Con_results/BayesNetInf/", file_dir, sep=""))
result_file_a = paste("C:/Users/Ravi/Desktop/Con_results/BayesNetInf/", file_dir, "/",result_file_a, sep="")

file_dir = "Feb62013"
file_dir = "Apr82013"
file_dir = "Apr182013"
result_file_a <- list.files(paste("C:/Users/rgoyal/Desktop/Simulation_Results/NetworkInference/", file_dir, sep=""))
result_file_a = paste("C:/Users/rgoyal/Desktop/Simulation_Results/NetworkInference/", file_dir, "/",result_file_a, sep="")


p_value_a_dmm = c()
mean_a_dmm = c()
truth_a_dmm = c()
mean_den_dmm = c()
truth_den_dmm = c()

view_graphs = TRUE

for (i in c(c(1:5),c(16:20), c(31:35))) {
  for (i in c(1:9)) {
    
    load(result_file_a[i])
    
    #genetic_seq_diagnostic(Initial_Data[[2]],T)
    #uncategorized_diagnostic(P_a, G_a, G, P, Ia, Il, R, Initial_Data) #Diagnostics that do not fall in a specific category
    #accuracy_G_P_diagnostic(P_a, G_a, Initial_Data)
    
    assort_info = network_properties_diagnostic(G_a, burnin = 000, view_graphs = view_graphs)  #Assumes G_a[[1]] is the original true network
    p_value_a_dmm = c(p_value_a_dmm, assort_info[[1]])
    mean_a_dmm = c(mean_a_dmm, assort_info[[2]])
    truth_a_dmm = c(truth_a_dmm, assort_info[[3]])
    mean_den_dmm = c(mean_den_dmm, assort_info[[4]])
    truth_den_dmm = c(truth_den_dmm, assort_info[[5]])
  }
  
  plot(truth_a_dmm, mean_a_dmm, xlim = c(-.3,.3), ylim = c(-.3,.3),col=c(1,1,1,2,1,1,3,3,2,1,1,2,2,2,2),
       cex = 2, pch = 16,
       xlab = "True Degree Assortativity Coefficient in Network",
       ylab = "Posterior Mean of Degree Assortativity Coefficient",
       main = "Compare True Degree Assortativity Coefficient vs Posterior"
  )
  lines(lowess(truth_a_dmm, mean_a_dmm), col = "blue", lwd = 5)
  abline(0,1, lwd = 5)
  
  
  plot(truth_den_dmm, mean_den_dmm, xlim = c(50,200), ylim = c(50,200),col=c(rep(1,5),rep(2,5),rep(3,5)),
       plot(truth_den_dmm, mean_den_dmm, xlim = c(1000,1100), ylim = c(1000,1100),col=c(rep(1,5),rep(2,5),rep(3,5)),
            
            cex = 2, pch = 16,
            xlab = "True Number of Edges in Network",
            ylab = "Posterior Mean of Number of Edges",
            main = "Compare True Number of Edges vs Posterior"
       )
       legend('topleft', c("Degree Mixing - Prior",
                           "",
                           "Disassortative",
                           "Assortative",
                           "Random"), 
              col = c("white", "white",c(1:6)),
              text.col = "black", pch = c(16)
       )
       lines(lowess(truth_den_dmm, mean_den_dmm), col = "blue", lwd = 5)
       abline(0,1, lwd = 5)
       
       
       p_value_a = c()
       mean_a = c()
       truth_a = c()
       mean_den= c()
       truth_den = c()
       
       for (i in c(c(6:15),c(21:30),c(36:45))) {
         load(result_file_a[i])
         
         #genetic_seq_diagnostic(Initial_Data[[2]],T)
         #uncategorized_diagnostic(P_a, G_a, G, P, Ia, Il, R, Initial_Data) #Diagnostics that do not fall in a specific category
         #accuracy_G_P_diagnostic(P_a, G_a, Initial_Data)
         
         assort_info = network_properties_diagnostic(G_a, burnin = 5000)  #Assumes G_a[[1]] is the original true network
         p_value_a = c(p_value_a, assort_info[[1]])
         mean_a = c(mean_a, assort_info[[2]])
         truth_a = c(truth_a, assort_info[[3]])
         mean_den= c(mean_den, assort_info[[4]])
         truth_den = c(truth_den, assort_info[[5]])
       }
       
       plot(truth_a, mean_a, xlim = c(-.4,.4), ylim = c(-.4,.4),col=c(c(2,2,2,2,1,1,1,2,1,1),c(2,2,1,2,1,2,1,1,2,1)+2, c(1,1,2,2,1,2,2,1,2,1)+4),
            cex = 2, pch = 16,
            xlab = "True Degree Assortativity Coefficient in Network",
            ylab = "Posterior Mean of Degree Assortativity Coefficient",
            main = "Compare True Degree Assortativity Coefficient vs Posterior"
       )
       legend('topleft', c("Degree Mixing - Prior",
                           "",
                           "Disassortative - Weak", "Disassortative - Strong",
                           "Assortative - Weak", "Assortative - Strong",
                           "Random - Weak", "Random - Strong"), 
              col = c("white", "white",c(1:6)),
              text.col = "black", pch = c(16)
       )
       lines(lowess(truth_a, mean_a), col = "blue", lwd = 5)
       abline(0,1, lwd = 5)
       
       plot(truth_den, mean_den, xlim = c(125,175), ylim = c(125,175),
            col=c(c(2,2,2,2,1,1,1,2,1,1),c(2,2,1,2,1,2,1,1,2,1)+2, c(1,1,2,2,1,2,2,1,2,1)+4),
            cex = 2, pch = 16,
            xlab = "True Number of Edges in Network",
            ylab = "Posterior Mean of Number of Edges",
            main = "Compare True Number of Edges vs Posterior"
       )
       legend('topleft', c("Degree Mixing - Prior",
                           "",
                           "Disassortative - Weak", "Disassortative - Strong",
                           "Assortative - Weak", "Assortative - Strong",
                           "Random - Weak", "Random - Strong"), 
              col = c("white", "white",c(1:6)),
              text.col = "black", pch = c(16)
       )
       
       lines(lowess(truth_den, mean_den), col = "blue", lwd = 5)
       abline(0,1, lwd = 5)
       
       
       
       par(mfrow = c(2,3))
       for (i in c(10,11,12,14,15)) { #Disassort and Weak
         load(result_file_a[i])
         
         assort_info = network_properties_diagnostic(G_a, burnin = 1000)  #Assumes G_a[[1]] is the original true network
         
       }
       
       par(mfrow = c(2,3))
       for (i in c(23,25,27,28,30)) { #Assortative and Weak
         load(result_file_a[i])
         
         assort_info = network_properties_diagnostic(G_a, burnin = 1000)  #Assumes G_a[[1]] is the original true network
         
       }
       
       par(mfrow = c(2,3))
       for (i in c(36,37,40,43,45)) { #Random and Weak
         load(result_file_a[i])
         
         assort_info = network_properties_diagnostic(G_a, burnin = 1000)  #Assumes G_a[[1]] is the original true network
         
       }
       
       #######Strong Prior
       
       par(mfrow = c(2,3))
       for (i in c(6,7,8,9,13)) { #Disassort and Weak
         load(result_file_a[i])
         
         assort_info = network_properties_diagnostic(G_a, burnin = 1000)  #Assumes G_a[[1]] is the original true network
         
       }
       
       par(mfrow = c(2,3))
       for (i in c(21,22,24,26,29)) { #Assortative and Weak
         load(result_file_a[i])
         
         assort_info = network_properties_diagnostic(G_a, burnin = 1000)  #Assumes G_a[[1]] is the original true network
         
       }
       
       par(mfrow = c(2,3))
       for (i in c(38,39,41,42,44)) { #Random and Weak
         load(result_file_a[i])
         
         assort_info = network_properties_diagnostic(G_a, burnin = 1000)  #Assumes G_a[[1]] is the original true network
         
       }
       
       c(c(2,2,2,2,1,1,1,2,1,1),c(2,2,1,2,1,2,1,1,2,1)+2, c(1,1,2,2,1,2,2,1,2,1)+4),
       
       

