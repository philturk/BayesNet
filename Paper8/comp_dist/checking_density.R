# --- combined two files into one (sampling_sim and data_restructuring)

# --- Initialize data and objects
set.seed(2021)
library(sna)
source("Network.R")

#Independent sub-samples of network to take for each sampling proportion
nsims <- 1000

# Sub-sampling proportions to simulate
samp_prop <- seq(from = .9, to = 0.2, by = -0.1)

#Establish array to store simulated data for each simulation at each sampling proportion.
sim_dist <- array(NA, dim = c(length(samp_prop), 23, nsims), 
                  dimnames = list(samp_prop, 1:23, NULL))

# --- Density of Largest Component in Full Network
# comp_membership = component.dist(GeneticNetwork)$cdist[1:23] #can i use this to create a table where density is entered or 0 if no comp of that size
comp_membership = component.dist(GeneticNetwork)$membership
comp_largest_id = which(comp_membership %>% table() == max(comp_membership %>% table()))
which(comp_membership == comp_largest_id ) %>% length()

comp_largest_ids = which(comp_membership == comp_largest_id ) 
comp_largest_net = get.inducedSubgraph(GeneticNetwork, v = comp_largest_ids)
summary(comp_largest_net ~ edges) #97
summary(comp_largest_net ~ density) #38.33% 

# ---Add Density Checking into Loop for different component sizes and 
for(i in 1:length(samp_prop)){
  sample_rate <- samp_prop[i]
  for(j in 1:nsims){
    #Randomly sample vertices to delete
    deleted_vertices <- sample(1:1258, round(1258*(1-sample_rate)), replace = FALSE)
    
    #Build subnetwork without deleted vertices
    sampledNetwork <- network(GeneticNetwork[-deleted_vertices, -deleted_vertices])
    
    # Compute cluster size distribution and store
    # sim_dist[i,,j] <- component.dist(sampledNetwork)$cdist[1:23]
    
    # compute densities instead of component distributions... at each size
  }
}

# ---Save the results
new_list <- plyr::adply(sim_dist,1,unlist,.id = NULL)
# write.csv(new_list, paste0(file_location, "//sim_dists-comp_dist-1000.csv"))
