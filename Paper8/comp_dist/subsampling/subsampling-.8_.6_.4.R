# --- FILE 1
set.seed(2021)

#Build complete genetic network
# it seems like the way that source() works is that it 
# adds together current working directory (can check using getwd())
# and whatever is in brackets here
getwd()
source("comp_dist//Network.R") 
# ^^ seems to only work when files are in subsample folder
# but we're referencing comp_dist folder?
# no idea what's going on with this 
# source("Network.R")

#Independent sub-samples of network to take for each sampling proportion
nsims <- 100

#Sub-sampling proportions to simulate
# sample proportion is number of ppl that were sampled
# samp_prop <- seq(from = .9, to = 0.2, by = -0.1)
# replace samp_prop that used (.9, .8 ... .2) with just .8, .6 and .4
# samp_prop <- c(.8, .6, .4)
# do for all three of these later
# for now just do for .8
samp_prop <- rep(.8, 100)

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
``
for(i in 1:length(samp_prop)){
  sample_rate <- samp_prop[i]
  for(j in 1:nsims){
    #Randomly sample vertices to delete
    deleted_vertices <- sample(1:1258, round(1258*(1-sample_rate)), replace = FALSE)
    
    #Build subnetwork without deleted vertices
    sampledNetwork <- network(GeneticNetwork[-deleted_vertices, -deleted_vertices])
    
    # Compute cluster size distribution and store
    sim_dist[i,,j] <- component.dist(sampledNetwork)$cdist[1:23]
  }
}
# created 814 columns of 23 - weird - then everything else is NA
test = sim_dist
colnames(test)
rownames(test) #these are our .8 repeated values
length(rownames(test)) #we get 100 ok
# so we want to replace with a diff set of rownames 8_1, 8_2...etc
list = c(1:100)
list2 = paste0("8_", list)
rownames(test) = list2
rownames(test)

# cool this worked and it's in the format where columns are simulations.
# and the first 23 rows are for .8, second 23 rows are for .8...etc
# new_list <- plyr::adply(sim_dist,1,unlist,.id = NULL)
new_list <- plyr::adply(test,1,unlist,.id = NULL)
file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
write.csv(new_list, paste0(file_location, "//sim_dists.8-comp_dist-100.csv"))


# --------------- GOAL OF DATA RESTRUCTURING ---------------
# BEFORE we want one df per proportion (.9, .8, ... .2) with 1000 sims
# BEFORE one line per count of components per simulation
# NOW we want one df per proportion (.8, .8, ... .8) with 100 sims <- then collapse them into one big DF with column 


# --------------- MELT ALL PROPORTIONS
# we've written a variety of these melt_all_props functions
# they have varied in the number of proportions (for .2, .3... .9 we have 8)
# and the number of simulations (originally we had 100 and then 1000)
# in our current case (for subsampling) we want 100 proportions (.8, .8..., .8)
# for 100 samples
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

# ----- CODE TO RUN
# load objects created in "sampling_sim.R"
sim_dist1 = read.csv(paste0(file_location, "//sim_dists.8-comp_dist-100.csv"))
df_all_sim = melt_all_props100by100(data = sim_dist1)#184000 rows = 1000*8*23, ok good

# check/look at it
df_all_sim[22998:23002,]
length(df_all_sim$X) #230000 = 100*100*23

# don't need to separate by proportion 
# so all 100 subsamples from each (of our same sized, .8, proportioned samples)
# is in the same file

# save each df into a separate csv
write.csv(df_all_sim, paste0(file_location, "//prop80//subsamples100by100.csv"), row.names = FALSE)
