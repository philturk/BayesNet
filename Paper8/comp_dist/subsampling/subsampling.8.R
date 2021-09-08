set.seed(2021)
library(sna)

#Build complete genetic network
getwd()
source("comp_dist//Network.R") 

nsims = 100 #number of subsamples for each sampling proportion
row_names = paste0("8_", c(1:100)) #row names for array
sample_rate = .8 #sampling rate
file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
size = 1258 # nodes(GeneticNetwork) # size of genetic network, 1258 in our case
samplesize = round(size*sample_rate) # size of sampled genetic network

melt_all_props100by100 = function(data){
  library(reshape2)
  # create one master df with the following columns:
  # X - the size of the component
  # variable - which sim it came from
  # value - count of components
  
  #refers to beg and end row lines for each prop of sim_dist
  beg_prop = seq(from = 1, to = 2300, by = 23) 
  end_prop = seq(from = 23, to = 2300, by = 23)
  # samp_prop <- rep(.8, 100) #does our function even use this?
  df = data.frame()
  df_all_sim = data.frame()
  for (i in 1:length(beg_prop)){
    beg = beg_prop[i]
    end = end_prop[i]
    df = as.data.frame(data[beg:end,])
    df$X = c(1:23) #leave as 1 to 23 because all props should have this component count
    df = melt(df, id.vars="X")
    df_all_sim = rbind(df_all_sim, df)
  }
  return(df_all_sim)
}

#Establish array to store simulated data for each simulation at each sampling proportion.
sim_dist <- array(NA, dim = c(length(row_names), 23, nsims), 
                  dimnames = list(row_names, 1:23, NULL))

# to try to keep the same networks as our previous, get deleted_vertices for all 
# for(i in 1:length(row_names)){
#   del <- sample(1:size, round(size*(1-sample_rate)), replace = FALSE)
#   deleted_vertices = append(deleted_vertices, del)
# }
# deleted_vertices[1]

for(i in 1:length(row_names)){
  deleted_vertices <- sample(1:size, round(size*(1-sample_rate)), replace = FALSE)
  
  sampledNetwork <- network(GeneticNetwork[-deleted_vertices, -deleted_vertices])
  
  for(j in 1:nsims){
    #Randomly sample vertices to delete
    deleted_vertices <- sample(1:samplesize, round(samplesize*(1-sample_rate)), replace = FALSE)
    
    #Build subnetwork without deleted vertices
    subsampledNetwork <- network(sampledNetwork[-deleted_vertices, -deleted_vertices])
    
    # Compute cluster size distribution and store
    sim_dist[i,,j] <- component.dist(subsampledNetwork)$cdist[1:23]
  }
  print(paste0("completed sim", as.character(i), " out of ", as.character(length(row_names))))
}

# change format and save 
new_list <- plyr::adply(sim_dist,1,unlist,.id = NULL)
# note to self: can go into file & sumproduct to check for correct # of nodes 
write.csv(new_list, paste0(file_location, "//sim_dists.8-comp_dist-100.csv"))

# change format and save
# new_list = read.csv(paste0(file_location, "//sim_dists.8-comp_dist-100.csv"))
df_all_sim = melt_all_props100by100(data = new_list)
write.csv(df_all_sim, paste0(file_location, "//prop80//subsamples100by100.csv"), row.names = FALSE)

# check/look at it
df_all_sim[22998:23002,]
length(df_all_sim$X) #230000 = 100*100*23