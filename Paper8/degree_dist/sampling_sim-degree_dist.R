library(sna)
library(janitor)
set.seed(2021)

#Build complete genetic network
file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(file_location, "//Network.R"))

#Independent sub-samples of network to take for each sampling proportion
nsims <- 100

#Sub-sampling proportions to simulate
# sample proportion is number of ppl that were sampled
samp_prop <- seq(from = .95, to = 0.2, by = -0.05)

#Compute cluster size distribution for complete network
# what we would have if we sampled everyone in the pop
true_dist <- unname(unlist(tabyl(factor(degree(GeneticNetwork, gmode="graph"), levels=0:20))[2]))

#Establish array to store simulated data for each simulation at each sampling proportion.
#Can assume degrees only get smaller, so max of size = 21 is sufficient (max size in full net)
# if we didn't know true dist and we picked diff sample sizes what
# would we infer the population looks like
sim_dist <- array(NA, dim = c(length(samp_prop), 21, nsims), 
                  dimnames = list(samp_prop, 1:21, NULL))

for(i in 1:length(samp_prop)){
  sample_rate <- samp_prop[i]
  for(j in 1:nsims){
    #Randomly sample vertices to delete
    deleted_vertices <- sample(1:1258, round(1258*(1-sample_rate)), replace = FALSE)
    
    #Build subnetwork without deleted vertices
    sampledNetwork <- network(GeneticNetwork[-deleted_vertices, -deleted_vertices])
    
    # put degree distribution in factor form to add zeros where count=0
    # dd = factor(degree(sampledNetwork, gmode="graph"), levels=0:20)
    # dd = tabyl(dd)[2]
    # dd = unname(unlist(dd))

    #Store degree distribution in object
    sim_dist[i,,j] <- unname(unlist(tabyl(factor(degree(sampledNetwork, gmode="graph"), levels=0:20))[2]))
  }
}
#Save the results
#save()

new_list <- plyr::adply(sim_dist,1,unlist,.id = NULL)
write.csv(new_list, paste0(file_location, "//sim_dists-degree_dist.csv"))
write.csv(true_dist, paste0(file_location, "//true_dist-degree_dist.csv"))