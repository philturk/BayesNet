args = (commandArgs(TRUE))
eval(parse(text = args[[1]]))
eval(parse(text = args[[2]]))
cat("num_samples = ", ns, "1st rep = ", seed)

library(parallel)
source("BIv3_densSampInvPropDeg_func.R")

#Number of cores to use
no_cores = 30

#Set random seeds as 1000 + simno
rand_seed = (1000 + seed):(1000 + seed + no_cores-1)
#Start cluster
cl <- makeCluster(no_cores)

#Run replicates of function on cluster
parLapply(cl, rand_seed, dens_biased_sim, num_samples = ns)

# Close cluster
stopCluster(cl)
