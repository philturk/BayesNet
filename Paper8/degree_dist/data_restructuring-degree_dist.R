library(reshape2)
# WARNING
# code has not been completely switched from comp_dist to degree_dist yet

file_location = dirname(rstudioapi::getActiveDocumentContext()$path)

# --------------- GOAL OF DATA RESTRUCTURING ---------------
# we want one df per proportion (.95, .9, ... .2)
# where each degree has it's own line
# and another column specifies which simulation it's from

# --------------- FUNCTIONS ---------------
# --------------SPLIT ONE SIMULATION
split_one_sim = function(beg = 1, end = 23, dataframe = NULL){
  # beg is starting row in df to split
  # end is last row in df to split
  # df is the dataframe we're splitting
  
  # initiate empty objects
  row_count = c()
  data = data.frame(simulation = character(), comp_size=integer())
  # loop to go over each degree in df
  for (i in beg:end){
    rows = dataframe$value[i] #how many rows to create
    input = dataframe$X[i] #which value to input into row
    if (rows > 0){
      start = sum(row_count)
      data[(start+1):(start+rows),2] = input
      row_count = append(row_count, rows)
    }
  }
  
  # add a column with which sim it's from
  var = paste0("sim_", (beg-1)/23+1)
  data[,1] = var
  
  return(data)
}

# --------------- COMBINE ALL SIMULATIONS
combine_all_sim = function(dataframe=NULL){
  # dataframe input should have 2300 lines and each 
  # 23 lines a new simulation results should start
  beg_100 = seq(from = 1, to = 2300, by = 23)
  end_100 = seq(from = 23, to = 2300, by = 23)
  combined = data.frame(simulation = character(), comp_size=integer())
  for (i in 1:length(beg_100)){
    beg = beg_100[i]
    end = end_100[i]
    data = split_one_sim(beg=beg, end=end, dataframe = dataframe)
    combined = rbind(combined, data)
  }
  return(combined)
}

# --------------- MELT ALL PROPORTIONS
melt_all_props = function(data = sim_dist_100){
  library(reshape2)
  # create one master df with the following columns:
  # X - the degree
  # variable - which sim it came from
  # value - count of degrees
  
  beg_prop = seq(from = 1, to = 346, by = 23)
  end_prop = seq(from = 23, to = 368, by = 23)
  samp_prop = seq(from = .95, to = 0.2, by = -0.05)
  df = data.frame()
  df_all_sim = data.frame()
  for (i in 1:length(beg_prop)){
    beg = beg_prop[i]
    end = end_prop[i]
    df = as.data.frame(data[beg:end,])
    df$X = c(1:23) #leave as 1 to 23 because all props should have this degree count
    df = melt(df, id.vars="X")
    df_all_sim = rbind(df_all_sim, df)
  }
  return(df_all_sim)
}

# ----- CODE TO RUN
# load objects created in "sampling_sim.R"
sim_dist_100 = read.csv(paste0(file_location, "//sim_dists-degree_dist.csv"))

# melt dataframe with all proportions and simulations into one master df_all_sim
df_all_sim = melt_all_props(data = sim_dist_100)
# define start and end rows in df_all_sim that separate each proportion
start_row = seq(from = 1, to = 36800, by = 2300)
end_row = seq(from = 2300, to = 36800, by = 2300)

data_95 = df_all_sim[start_row[1]:end_row[1],3]

# create a separate df for each proportion
df_95 = combine_all_sim(dataframe = df_all_sim[start_row[1]:end_row[1],])
df_9 = combine_all_sim(dataframe = df_all_sim[start_row[2]:end_row[2],]) 
df_85 = combine_all_sim(dataframe = df_all_sim[start_row[3]:end_row[3],])
df_8 = combine_all_sim(dataframe = df_all_sim[start_row[4]:end_row[4],]) 
df_75 = combine_all_sim(dataframe = df_all_sim[start_row[5]:end_row[5],])
df_7 = combine_all_sim(dataframe = df_all_sim[start_row[6]:end_row[6],]) 
df_65 = combine_all_sim(dataframe = df_all_sim[start_row[7]:end_row[7],])
df_6 = combine_all_sim(dataframe = df_all_sim[start_row[8]:end_row[8],]) 
df_55 = combine_all_sim(dataframe = df_all_sim[start_row[9]:end_row[9],])
df_5 = combine_all_sim(dataframe = df_all_sim[start_row[10]:end_row[10],]) 
df_45 = combine_all_sim(dataframe = df_all_sim[start_row[11]:end_row[11],])
df_4 = combine_all_sim(dataframe = df_all_sim[start_row[12]:end_row[12],]) 
df_35 = combine_all_sim(dataframe = df_all_sim[start_row[13]:end_row[13],])
df_3 = combine_all_sim(dataframe = df_all_sim[start_row[14]:end_row[14],]) 
df_25 = combine_all_sim(dataframe = df_all_sim[start_row[15]:end_row[15],])
df_2 = combine_all_sim(dataframe = df_all_sim[start_row[16]:end_row[16],])
df_only_95 = df_all_sim[start_row[1]:end_row[1],]

# save each df into a separate csv
write.csv(df_95, paste0(file_location, "//proportions//prop95-dd.csv"), row.names = FALSE)
write.csv(df_9, paste0(file_location, "//proportions//prop90-dd.csv"), row.names = FALSE)
write.csv(df_85, paste0(file_location, "//proportions//prop85-dd.csv"), row.names = FALSE)
write.csv(df_8, paste0(file_location, "//proportions//prop80-dd.csv"), row.names = FALSE)
write.csv(df_75, paste0(file_location, "//proportions//prop75-dd.csv"), row.names = FALSE)
write.csv(df_7, paste0(file_location, "//proportions//prop70-dd.csv"), row.names = FALSE)
write.csv(df_65, paste0(file_location, "//proportions//prop65-dd.csv"), row.names = FALSE)
write.csv(df_6, paste0(file_location, "//proportions//prop60-dd.csv"), row.names = FALSE)
write.csv(df_55, paste0(file_location, "//proportions//prop55-dd.csv"), row.names = FALSE)
write.csv(df_5, paste0(file_location, "//proportions//prop50-dd.csv"), row.names = FALSE)
write.csv(df_45, paste0(file_location, "//proportions//prop45-dd.csv"), row.names = FALSE)
write.csv(df_4, paste0(file_location, "//proportions//prop40-dd.csv"), row.names = FALSE)
write.csv(df_35, paste0(file_location, "//proportions//prop35-dd.csv"), row.names = FALSE)
write.csv(df_3, paste0(file_location, "//proportions//prop30-dd.csv"), row.names = FALSE)
write.csv(df_25, paste0(file_location, "//proportions//prop25-dd.csv"), row.names = FALSE)
write.csv(df_2, paste0(file_location, "//proportions//prop20-dd.csv"), row.names = FALSE)
write.csv(df_all_sim, paste0(file_location, "//proportions//df_all_sim-dd.csv"), row.names = FALSE)
write.csv(df_only_95, paste0(file_location, "//proportions//df_only_95-dd.csv"), row.names = FALSE)
# call to map function (similar to apply, takes a function & vector & return list/vector)
# partial for map - sets default inputs

# # --- check
# # for check
# samp_prop <- seq(from = .95, to = 0.2, by = -0.05)
# for(i in 1:length(samp_prop)){
#   sample_rate <- samp_prop[i]
#   print(round(1258*(sample_rate)))
#   }
# # 1195, 1132, 1069, 1006, 944, 881, 818, 755, 692, 629, 566, 503, 440, 377, 314, 252
#   
# # check the sim is sampling .95 nodes versus components
# crossprod(true_dist_1$X,true_dist_1$x)*100 #125800
# crossprod(df_all_sim$X[1:2300],df_all_sim$value[1:2300]) #119500
# 125800*.95 #119510 <-- shouldn't this be the crossproduct of above? off by 10
# # makes sense because this is with 100 sims 
# # but without we needed it to be the nearest node
# 
# crossprod(true_dist_1$X,true_dist_1$x)*100 #125800
# crossprod(df_all_sim$X[2301:4600],df_all_sim$value[2301:4600]) #113200
# 125800*.9 #113220 <-- off by 10 again
# resolved! got rounded. see code above

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

# --- checks for output of combine_all__sim
# sum(df_95$comp_size) #same node count
# crossprod(df_all_sim[start_row[1]:end_row[1],]$X, df_all_sim[start_row[1]:end_row[1],]$value)
# 
# length(df_95$comp_size)
# sum(df_all_sim[start_row[1]:end_row[1],]$value)
# 
# # --- checks for output of combine_all__sim
# sum(df_95$comp_size) #same node count
# crossprod(df_all_sim[start_row[1]:end_row[1],]$X, df_all_sim[start_row[1]:end_row[1],]$value)
# 
# length(df_95$comp_size)
# sum(df_all_sim[start_row[1]:end_row[1],]$value)
# 
# sum(df_2$comp_size) #same node count
# crossprod(df_all_sim[start_row[16]:end_row[16],]$X, df_all_sim[start_row[16]:end_row[16],]$value)
# 
# length(df_2$comp_size)
# sum(df_all_sim[start_row[16]:end_row[16],]$value)
# 
# # other checks
# head(df_all_sim[start_row[1]:end_row[1],], 5) #ok they're def different
# head(df_all_sim[start_row[2]:end_row[2],], 5) #ok they're def different
# start_row[1]
# end_row[1]
# start_row[2]
# end_row[2]
# 
# tail(df_all_sim[start_row[2]:end_row[2],], 10) #this does not add up?
# tail(df_9, 10) #shows 100th sim has 21 * 2 and 18 * 1, but above line doesn't
# 
# sum(df_all_sim[start_row[2]:end_row[2],]$value) #76702  <- this should be the # of observations
# 
