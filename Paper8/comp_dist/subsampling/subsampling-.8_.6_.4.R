# --- combined two files into one (sampling_sim and data_restructuring)

# --- FILE 1
set.seed(2021)

#Build complete genetic network
file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(file_location, "//Network.R"))
source("Network.R")

#Independent sub-samples of network to take for each sampling proportion
nsims <- 1000

#Sub-sampling proportions to simulate
# sample proportion is number of ppl that were sampled
samp_prop <- seq(from = .9, to = 0.2, by = -0.1)

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
    
    # Compute cluster size distribution and store
    sim_dist[i,,j] <- component.dist(sampledNetwork)$cdist[1:23]
  }
}
# created 814 columns of 23 - weird - then everything else is NA

#Save the results
#save()
# cool this worked and it's in the format where columns are simulations.
# and the first 23 rows are for .9, second 23 rows are for .80...etc
new_list <- plyr::adply(sim_dist,1,unlist,.id = NULL)
write.csv(new_list, paste0(file_location, "//sim_dists-comp_dist-1000.csv"))


# --- FILE 2

file_location = dirname(rstudioapi::getActiveDocumentContext()$path)

# --------------- GOAL OF DATA RESTRUCTURING ---------------
# we want one df per proportion (.9, .8, ... .2) with 1000 sims
# one line per count of components per simulation

# --------------- MELT ALL PROPORTIONS
melt_all_props1000 = function(data = sim_dist_100){
  library(reshape2)
  # create one master df with the following columns:
  # X - the size of the component
  # variable - which sim it came from
  # value - count of components
  
  # beg_prop = seq(from = 1, to = 346, by = 23) #(346-1)/23 = 15 = (16-1)
  # end_prop = seq(from = 23, to = 368, by = 23) #368/23=16 props but now we have 8 (half of that)
  beg_prop = seq(from = 1, to = 162, by = 23) #just guessing 23*(8-1) = 161
  end_prop = seq(from = 23, to = 184, by = 23) #changed it to 8*23=184 - double check
  samp_prop = seq(from = .9, to = 0.2, by = -0.1)
  df = data.frame()
  df_all_sim = data.frame()
  for (i in 1:length(beg_prop)){ #length is 8
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
sim_dist_1000 = read.csv(paste0(file_location, "//sim_dists-comp_dist-1000.csv"))
max_comp = 23
# melt dataframe with all proportions and simulations into one master df_all_sim
df_all_sim = melt_all_props1000(data = sim_dist_1000)#184000 rows = 1000*8*23, ok good

df_all_sim[22998:23002,]

# define start and end rows in df_all_sim that separate each proportion
start_row = seq(from = 1, to = length(df_all_sim$X), by = nsims*max_comp)
end_row = seq(from = nsims*max_comp, to = length(df_all_sim$X), by = nsims*max_comp)

# create a separate df for each proportion
df_9 = df_all_sim[start_row[1]:end_row[1],]
df_8 = df_all_sim[start_row[2]:end_row[2],] 
df_7 = df_all_sim[start_row[3]:end_row[3],]
df_6 = df_all_sim[start_row[4]:end_row[4],] 
df_5 = df_all_sim[start_row[5]:end_row[5],]
df_4 = df_all_sim[start_row[6]:end_row[6],] 
df_3 = df_all_sim[start_row[7]:end_row[7],]
df_2 = df_all_sim[start_row[8]:end_row[8],]

# save each df into a separate csv
write.csv(df_9, paste0(file_location, "//props1000//df_90-comp_dist-1000.csv"), row.names = FALSE)
write.csv(df_8, paste0(file_location, "//props1000//df_80-comp_dist-1000.csv"), row.names = FALSE)
write.csv(df_7, paste0(file_location, "//props1000//df_70-comp_dist-1000.csv"), row.names = FALSE)
write.csv(df_6, paste0(file_location, "//props1000//df_60-comp_dist-1000.csv"), row.names = FALSE)
write.csv(df_5, paste0(file_location, "//props1000//df_50-comp_dist-1000.csv"), row.names = FALSE)
write.csv(df_4, paste0(file_location, "//props1000//df_40-comp_dist-1000.csv"), row.names = FALSE)
write.csv(df_3, paste0(file_location, "//props1000//df_30-comp_dist-1000.csv"), row.names = FALSE)
write.csv(df_2, paste0(file_location, "//props1000//df_20-comp_dist-1000.csv"), row.names = FALSE)