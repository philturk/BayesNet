library(reshape2)

# we want one (mixed) model per % of the original sim
# one model for .95, .9.... etc (so one df for each)
# so we want a df where each component has it's own line in the df

# load objects created in "sampling_sim.R"
sim_dist_100 = read.csv("C://Code//msu//gra//BayesNet//Paper8//sim_dists.csv")
true_dist_1 = read.csv("C://Code//msu//gra//BayesNet//Paper8//true_dist.csv")

# keep until discrepancy of 10 nodes in loop is solved
# delete after
df = as.data.frame(sim_dist_100[1:23,])
df = melt(df, id.vars="X")
head(df)

# --------------- FUNCTIONS

# --------------SPLIT ONE SIMULATION
split_one_sim = function(beg = 1, end = 23, dataframe = df){
  # beg is starting row in df to split
  # end is last row in df to split
  # df is the dataframe we're splitting
  
  # initiate empty objects
  row_count = c()
  data = data.frame(simulation = character(), comp_size=integer())
  # loop to go over each component size in df
  for (i in beg:end){
    rows = dataframe$value[i] #how many rows to create
    input = dataframe$X[i] #which value to input into row
    if (rows > 0){
      start = sum(row_count)
      data[(start+1):(start+rows),2] = input #I think if there's a prob it's in this line
      row_count = append(row_count, rows)
    }
  }
  
  # add a column with which sim it's from
  var = paste0("sim_", (beg-1)/23+1)
  data[,1] = var
  
  return(data)
}

# --------------- COMBINE ALL SIMULATIONS
combine_all_sim = function(){
  beg_100 = seq(from = 1, to = 2300, by = 23)
  end_100 = seq(from = 23, to = 2300, by = 23)
  combined = data.frame(simulation = character(), comp_size=integer())
  for (i in 1:length(beg_100)){
    beg = beg_100[i]
    end = end_100[i]
    data = split_one_sim(beg=beg, end=end)
    combined = rbind(combined, data)
  }
  return(combined)
}

# --------------- MELT ALL PROPORTIONS
melt_all_props = function(){
  beg_prop = seq(from = 1, to = 346, by = 23)
  end_prop = seq(from = 23, to = 368, by = 23)
  samp_prop = seq(from = .95, to = 0.2, by = -0.05)
  # start_row = seq(from = 0, to = 34500, by = 2300) #maybe don't need if I use rbind
  df = data.frame()
  df_all_sim = data.frame()
  for (i in 1:length(beg_prop)){
    beg = beg_prop[i]
    end = end_prop[i]
    df = as.data.frame(sim_dist_100[beg:end,])
    df$X = c(1:23) #leave as 1 to 23 because all props should have this component count
    df = melt(df, id.vars="X")
    df_all_sim = rbind(df_all_sim, df)
  }
  return(df_all_sim)
}


# last step - nest them all in each other
df_all_sim = melt_all_props()




# --- check
# WHY ISN'T THIS ADDING UP?
# check the sim is sampling .95 nodes versus components
crossprod(true_dist_1$X,true_dist_1$x)*100 #125800
crossprod(df_all_sim$X[1:2300],df_all_sim$value[1:2300]) #119500
125800*.95 #119510 <-- shouldn't this be the crossproduct of above? off by 10

crossprod(true_dist_1$X,true_dist_1$x)*100 #125800
crossprod(df_all_sim$X[2301:4600],df_all_sim$value[2301:4600]) #113200
125800*.9 #113220 <-- off by 10 again




# can build in print-outs for checks if desired...
# # ---checks for output of one loop of split_one_sim()
# test_1 = split_one_sim(beg=1, end=23)
# df$value[1:23]
# table(test_1$comp_size)
# 
# sum(df$value[1:23]) #798 components
# length(test_1$comp_size)#798 components YES
# 
# crossprod(df$X[1:23],df$value[1:23]) #1195 total nodes
# sum(test_1$comp_size) #1195 total nodes

# --- checks for output of combine_all_sim() for one proportion
# data = combine_all_sim() #length of data is 79938
# crossprod(df$X,df$value) #119500 total nodes
# sum(data$comp_size) #119500

