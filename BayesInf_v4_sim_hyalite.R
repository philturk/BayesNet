#WRAPPER SCRIPT TO SIMULATE ONE EPIDEMIC & PERFORM BAYESIAN INFERENCE - ON HYALITE CLUSTER
#create a Simulations directory before running (below)

####################Setup####################

#start time
start <- Sys.time()

#pass task ID from job array to R
args <- commandArgs(trailingOnly=TRUE)
id <- as.numeric(args[[1]])

####################Simulation####################

#call simulation script
source('BayesInf_v4_sim.R')

####################Results####################

#save seed, if simulation fails (no epidemic)
save(init_seed, file=paste("Simulations/simulation-",formatC(id, width=3, format="d", flag="0"),".rda", sep=""))

#create output object
output <- list(NULL)

#add simulation results to output
output[[3]] <- results_list

#add ID number to output
output[[1]] <- id

#end time
end <- Sys.time()

#run time
runtime <- end-start

#add runtime to output
output[[2]] <- runtime

#save simulation output with ID number (saves over seed)
save(output, file=paste("Simulations/simulation-",formatC(id, width=3, format="d", flag="0"),".rda", sep=""))

#output (list of 3):
#[[1]] simulation ID
#[[2]] run time
#[[3]] results
