#Build complete genetic network
source("Network.R")

#Independent sub-samples of network to take for each sampling proportion
nsims <- 100

#Sub-sampling proportions to simulate
# sample proportion is number of ppl that were sampled
samp_prop <- seq(from = .95, to = 0.2, by = -0.05)

#Compute cluster size distribution for complete network
library(sna)
# what we would have if we sampled everyone in the pop
true_dist <- component.dist(GeneticNetwork)$cdist[1:23]

#Establish array to store simulated data for each simulation at each sampling proportion.
#Can assume clusters only get smaller, so max of size = 23 is sufficient (max size in full net)
# if we didn't know true dist and we picked diff sample sizes what
# would we infer the population looks like
sim_dist <- array(NA, dim = c(length(samp_prop), 23, nsims), 
                  dimnames = list(samp_prop, 1:23, NULL))

for(i in 1:length(samp_prop)){
  sample_rate <- samp_prop[i]
  for(j in 1:nsims){
    #Randomly sample vertices to delete
    deleted_vertices <- sample(1:1258, round(1258*(1-sample_rate)), replace = FALSE)
    
    #Build subnetwork without deleted vertices
    sampledNetwork <- network(GeneticNetwork[-deleted_vertices, -deleted_vertices])
    
    #Compute cluster size distribution and store
    sim_dist[i,,j] <- component.dist(sampledNetwork)$cdist[1:23]
  }
}

#Save the results
#save()
# cool this worked and it's in the format where columns are simulations.
# and the first 23 rows are for .95, 
# second 23 rows are for .80...etc (20 1+368 rows total (includes header))
new_list <- plyr::adply(sim_dist,1,unlist,.id = NULL)
write.csv(new_list, "sim_dists.csv")
write.csv(true_dist, "true_dist.csv")