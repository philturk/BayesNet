set.seed(2021)
source("C://Code//msu//gra//BayesNet//Paper8//ergm_sims//gof_comp_size//ComponentGOF.R")

# note to self
# the input is a network object 
# use loop from previous work & save network object instead
# then loop over the saved networks and apply models below
# save outputs from models (gof plots, etc)


##################
# models that converged for loop over network objects
##################
model.01 = ergm(GeneticNetwork~edges)
model.02 = ergm(GeneticNetwork~edges + degree(1))
model.03 = ergm(GeneticNetwork~edges + degree(2)) #11 iterations
model.04 = ergm(GeneticNetwork~edges + degree(1) + degree(2)) 
model.05 = ergm(GeneticNetwork ~ gwdegree(decay=100,fixed=TRUE))

# loop inside loop to save relevant info about each model
summary = summary(model)
comp_fit = ComponentGOF(ergm=model, 10)
dd_fit = plot(gof(model))
diagnostics = mcmc.diagnostics(model) 

# after can do a visual of how coefficients change in the models as we go through
# want multiple samples per sampling rate?