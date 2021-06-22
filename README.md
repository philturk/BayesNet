# BayesNet

# There are 2 main files to execute the method (only select one):

## BayesInf_v4.R - Used for a single simulation run
## BayesInf_v4_sim.R - Used for the cluster

# The main files call the following programs:

## BayesInf_v4_func_read_params.R 
      - Contains all of the specificians for the simulations.
      - States type of network property of interest (e.g., mixing, degree distribution)
      
## BayesInf_v4_func_gen_epi_data.R 
      - Generates the ground truth 
      - Generates simulated datasets (genetic, behavioral survey, clinical data)
      
## BayesInf_v4_func_gen_genetic_data.R
      - A sub script for BayesInf_v4_func_gen_epi_data.R to simulate genetic data
      
## BayesInf_v4_func_initial_bayes_inf.R
      - The Bayesian inference needs an initial contact network and transmission network
      - Makes the initial transmission network consistent with the observed clinical data
      
## BayesInf_v4_func_bayes_inf.R
      - Updates the contact network, transmission network, and parameter estimates
      - Updating the contact network uses the CCMnet package (which uses python)
      
## BayesInf_v4_func_diagnostics.R
      - Mostly junk for now
      
      
# Relationship between the files:
      - BayesInf_v4.R / BayesInf_v4_sim.R calls all the files 
      - First BayesInf_v4_func_read_params.R
      - Second BayesInf_v4_func_gen_epi_data.R
          - BayesInf_v4_func_gen_epi_data.R calls BayesInf_v4_func_gen_genetic_data.R
      - Third BayesInf_v4_func_initial_bayes_inf.R
      - Fourth BayesInf_v4_func_bayes_inf.R
      
      

