set.seed(2021)
library(sna)

#Build complete genetic network
getwd()
source("comp_dist//Network.R") 

# declare variables
nsims = 100 #number of subsamples for each sampling proportion
row_names = paste0("6_", c(1:100)) #row names for array
sample_rate = .6 #sampling rate
file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
sim_dist <- array(NA, dim = c(length(row_names), 23, nsims), 
                  dimnames = list(row_names, 1:23, NULL)) #array to store sims

# function for later
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


for(i in 1:length(row_names)){
  for(j in 1:nsims){
    #Randomly sample vertices to delete
    deleted_vertices <- sample(1:1258, round(1258*(1-sample_rate)), replace = FALSE)
    
    #Build subnetwork without deleted vertices
    sampledNetwork <- network(GeneticNetwork[-deleted_vertices, -deleted_vertices])
    
    # Compute cluster size distribution and store
    sim_dist[i,,j] <- component.dist(sampledNetwork)$cdist[1:23]
  }
  print(paste0("completed sim", as.character(i), " out of ", as.character(length(row_names))))
}

# change format, save
new_list <- plyr::adply(sim_dist,1,unlist,.id = NULL)
write.csv(new_list, paste0(file_location, "//sim_dists.6-comp_dist-100.csv"))

# change format again, save
df_all_sim = melt_all_props100by100(data = new_list)
write.csv(df_all_sim, paste0(file_location, "//prop60//subsamples100by100.csv"), row.names = FALSE)

# check/look at it
df_all_sim[22998:23002,]
length(df_all_sim$X) #230000 = 100*100*23