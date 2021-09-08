set.seed(2021)
library(sna)

#Build complete genetic network
getwd()
source("comp_dist//Network.R") 

# description I wrote down: take one sample at 50%
# subsample 100 times at each of .45, .4... .2
# fit exp 3 param
# repeat 100 times

# declare variables to go from 50 to 20 in increments of .05 with 100 subsamples at each level
nsims = 100 #number of subsamples for each sampling proportion
row_names = seq(from=.45, to=.2, by=-.05) #row names for array
sample_rate = .5
file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
size = 1258 # nodes(GeneticNetwork) # size of genetic network, 1258 in our case
sim_dist <- array(NA, dim = c(length(row_names), 23, nsims), 
                  dimnames = list(row_names, 1:23, NULL)) #array to store sims


for(i in 1:length(row_names)){
  
  deleted_vertices <- sample(1:size, round(size*(1-sample_rate)), replace = FALSE)
  
  sampledNetwork <- network(GeneticNetwork[-deleted_vertices, -deleted_vertices])
  
  samplesize = round(size*sample_rate) # number of nodes for that sample_rate for subsample
  # will be the same for all 100*6 subsample rates times (since we repeat .5 for)
  
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

# change format, save
new_list <- plyr::adply(sim_dist,1,unlist,.id = NULL)
write.csv(new_list, paste0(file_location, "//sim_dists-incre-subsamp-comp_dist-100.csv"))

# function for later
melt_all_props7by100 = function(data){
  library(reshape2)
  # create one master df with the following columns:
  # X - the size of the component
  # variable - which sim it came from
  # value - count of components
  
  #refers to beg and end row lines for each prop of sim_dist
  beg_prop = seq(from = 1, to = 161, by = 23) 
  end_prop = seq(from = 23, to = 161, by = 23)
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


# change format again, save
df_all_sim = melt_all_props7by100(data = new_list)
write.csv(df_all_sim, paste0(file_location, "//incremental//subsamples70by100.csv"), row.names = FALSE)

# check/look at it
length(df_all_sim$X) #16 100 = 100 (subsamples) *7 (props) *23 (comp sizes)