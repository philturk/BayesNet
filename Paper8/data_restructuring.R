library(reshape2)

# we want one (mixed) model per % of the original sim
# one model for .95, .9.... etc (so one df for each)
# so we want a df where each component has it's own line in the df

# load objects created in "sampling_sim.R"
sim_dist_100 = read.csv("C://Code//msu//gra//BayesNet//Paper8//sim_dists.csv")
true_dist_1 = read.csv("C://Code//msu//gra//BayesNet//Paper8//true_dist.csv")

# --- first proportion. 
# After this will be in a loop where rows 1:23 are .95, 24:46 are .9
# see below for attempted loop
# create df in format: 
# one column for component size count, 
# one column for which sim it came from
df = as.data.frame(sim_dist_100[1:23,])
df = melt(df, id.vars="X")
head(df)

# this is a loop for each of the proportions... 
# create df with one column per proportion (includes 100 simulations)
for (i in c(1:16)){
  df = as.data.frame(sim_dist_100[(i+((i-1)*22)):(i*23),])
  df$X = c(1:23)
  df = melt(df, id.vars="X")
  df_all_sim[as.character(samp_prop[i])] = df$value
}

# --- split one sim


# # ------------create df with one row = 1 component
# # create a dataframe for one simulation (ex: X1)
# # will need to nest this in another loop to do the same for all of X1, X2... X100
# # delete after
# row_count = c()
# data = data.frame(comp_size=integer())
# for (i in 1:23){
#   rows = df$value[i] #how many rows to create
#   input = df$X[i] #which value to input into row
#   start = sum(row_count)
#   data[(start+1):(start+rows),1] = input
#   row_count = append(row_count, rows)
# }
# table(data$comp_size)
# ------------



# one sim from a df with one proportion
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

# # ---checks
# # AHHH why aren't these lining up?! -- stop to eat
# # create and combine all the sims from one proportion
# test_1 = split_one_sim(beg=1, end=23)
# df$value[1:23]
# table(test_1$comp_size)
# 
# sum(df$value[1:23]) #798 components
# length(test_1$comp_size)#798 components YES
# 
# crossprod(df$X[1:23],df$value[1:23]) #1195 total nodes
# sum(test_1$comp_size) #1195 total nodes
# 
# test2 = split_one_sim(beg=24, end=46)
# sum(df$value[24:46]) #802
# table(test2$comp_size)
# df$value[24:46] 
# 
# test3 = split_one_sim(beg=47, end=69)
# table(test3$comp_size)
# df$value[47:69]

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

data = combine_all_sim() #length of data is 79938
crossprod(df$X,df$value) #119500 total nodes
sum(data$comp_size) #119500
# 
# 
# sum(df$value) #79868
# table(data$comp_size[1:799])
# table(data$comp_size[800:800+802])
# beg_100
# end_100
# # 79938-sum(df$value) is 70
# 
