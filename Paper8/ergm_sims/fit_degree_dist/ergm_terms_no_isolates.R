# do the same thing the ergm_terms.R file is doing but with isolated nodes dropped. 
# Can it fit the comp distribution better without this extreme value?

set.seed(2021)
# table of common ergm terms:
# https://statnet.org/nme/d2-ergmterms.html

# source("C://Code//msu//gra//BayesNet//Paper8//degree_dist//Network.R")
# source("Network.R")
source("C://Code//msu//gra//BayesNet//Paper8//ergm_sims//gof_comp_size//ComponentGOF.R")
# source("ComponentGOF.R")

# test = multiplex::rm.isol(GeneticNetwork)
set.vertex.attribute(GeneticNetwork, "Degree", sna::degree(GeneticNetwork, gmode = "graph"))
test = get.inducedSubgraph(GeneticNetwork, which(GeneticNetwork %v% "Degree" > 0))

#################################################
# edges + isolatededges + concurrent() + isolates
#################################################
ergm(test ~ edges + isolatededges + concurrent() + isolates)
# degenerate
# In addition: Warning message:
# In ergm_MCMC_sample(s, control, theta = mcmc.init, verbose = max(verbose -  :
# Unable to reach target effective size in iterations alotted.

#################################################
# edges + isolatededges + concurrent()
#################################################
ergm(test ~ edges + isolatededges + concurrent())
# degenerate

#################################################
# edges + isolatededges + concurrent()
#################################################
ergm(test ~ edges + concurrent())
# Model statistics 'concurrent' are not varying. This may indicate that the observed data occupies an extreme point in the sample space or that the estimation has reached a dead-end configuration.
# Optimizing with step length 0.0001.
# The log-likelihood improved by 0.0397.
# Estimating equations are not within tolerance region.
# Error in ergm.MCMLE(init, nw, model, initialfit = (initialfit <- NULL),  : 
# MCMLE estimation stuck. There may be excessive correlation between model terms, suggesting a poor model for the observed data. If target.stats are specified, try increasing SAN parameters.

#################################################
# edges
#################################################
model2.01 = ergm(test ~ edges)
summary(model2.01)
plot(gof(model2.01))
ComponentGOF(ergm=model2.01, 100)
mcmc.diagnostics(model2.01)


