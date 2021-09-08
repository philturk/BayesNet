set.seed(2021)
# table of common ergm terms:
# https://statnet.org/nme/d2-ergmterms.html

# source("C://Code//msu//gra//BayesNet//Paper8//degree_dist//Network.R")
# source("Network.R")
source("C://Code//msu//gra//BayesNet//Paper8//ergm_sims//gof_comp_size//ComponentGOF.R")

summary(GeneticNetwork)
?"ergm-terms"

# note to self
# To examine model diagnostics and check
# for degeneracy, use the mcmc.diagnostics() function.


##################
# edges only
##################
summary(GeneticNetwork ~ edges) #1020 edges in network
model.01 = ergm(GeneticNetwork~edges)
summary(model.01)
ComponentGOF(ergm=model.01, 10)
plot(gof(model.01))
mcmc.diagnostics(model.01) 
#none for mcmc.... why? because it's countable so don't need to estimate denominator?
# sizes obs   mean
# 1 659 253.56
# 2 102  38.67
# 3  24  14.22
# 4  15   5.67
# 5   9   3.22
# 6   4   1.67
# 7   3   0.78
# 8   3   0.67
# 9   4   0.33
# 10   0   0.33
# 11   1   0.11
# 12   2   0.00
# 13   0   0.11
# 14   1   0.11
# 15   0   0.00
# 16   0   0.11
# 17   0   0.00
# 18   0   0.00
# 19   0   0.00
# 20   1   0.00
# 21   1   0.00
# 22   0   0.00
# 23   1   0.00
# 24   0   0.00
# 25+   0   4.00


##################
# triangles only
##################
summary(GeneticNetwork ~ triangles) #1556 triangles in network
ergm(GeneticNetwork ~ triangles)
# error

##################
# triangles and edges
##################
ergm(GeneticNetwork ~ edges + triangles)
# error

##################
# degree(1)
##################
summary(GeneticNetwork ~ degree(1)) #237
ergm(GeneticNetwork~degree(1))
# error


##################
# degree(2)
##################
summary(GeneticNetwork ~ degree(2)) #100
ergm(GeneticNetwork~degree(2))
# error


##################
# edges + degree(1)
##################
model.02 = ergm(GeneticNetwork~edges + degree(1))
summary(model.02)
plot(gof(model.02))
ComponentGOF(ergm=model.02, 10)
mcmc.diagnostics(model.02)


##################
# edges + degree(2)
##################
model.03 = ergm(GeneticNetwork~edges + degree(2)) #11 iterations
summary(model.03)
plot(gof(model.03))
ComponentGOF(ergm=model.03, 10)
mcmc.diagnostics(model.03)



##################
# edges + degree(1) + degree(2)
##################
model.04 = ergm(GeneticNetwork~edges + degree(1) + degree(2)) 
#15 iterations first time, then 17 the second... what about seed?
summary(model.04)
mcmc.diagnostics(model.04)
plot(gof(model.04))
ComponentGOF(ergm=model.04, 10)
mcmc.diagnostics(model.04)


##################
# edges + degree(1) + degree(2) + degree(3)
##################
summary(GeneticNetwork~edges + degree(1) + degree(2) + degree(3))
model.06 = ergm(GeneticNetwork~edges + degree(1) + degree(2) + degree(3)) 
summary(model.06)
mcmc.diagnostics(model.06)
plot(gof(model.06))
ComponentGOF(ergm=model.06, 10)


##################
# gwdegree(decay=100,fixed=TRUE)
##################
summary(GeneticNetwork ~ gwdegree(decay=100,fixed=TRUE))
model.05 = ergm(GeneticNetwork ~ gwdegree(decay=100,fixed=TRUE))
plot(gof(model.05))
ComponentGOF(ergm=model.05, 10)
mcmc.diagnostics(model.05)



##################
# Kara's
# edges + gwdegree(decay=100,fixed=TRUE)
##################
summary(GeneticNetwork ~ edges + gwdegree(decay=100,fixed=TRUE)) #1020, 2040
ergm(GeneticNetwork~edges + gwdegree(decay=100,fixed=TRUE)) 
# this model takes FOREVER <-- abandoned ship


##################
# stars
##################
summary(GeneticNetwork ~ edges + kstar(1)) #1020, 2040
ergm(GeneticNetwork~edges + kstar(1)) 
# this model takes forever <-- abandoned ship
# notice: in undirected networks it says edges = kstar(1) so i'm repeating ops


##################
# isolates only
##################
summary(GeneticNetwork ~ isolates) #659
ergm(GeneticNetwork~isolates) 
# degenerate
# also got degenerate for 
# ergm(GeneticNetwork~edges + degree(0)) 
# which makes sense b/c it's the same thing


##################
# isolates & edges
##################
summary(GeneticNetwork ~ edges + isolates) #1020, 659
ergm(GeneticNetwork~edges + isolates) 
# degenerate


##################
# gwesp() only
##################
summary(GeneticNetwork ~ gwesp()) #659
ergm(GeneticNetwork~gwesp()) 
# degenerate

##################
# isolates & edges
##################
summary(GeneticNetwork ~ edges + gwesp()) #1020, 659
ergm(GeneticNetwork~edges + gwesp()) 
# degenerate

##################
# edges + concurrent
##################
model.06 = ergm(GeneticNetwork~ edges + concurrent())
summary(model.06)
plot(gof(model.06))
ComponentGOF(ergm=model.06, 100)
mcmc.diagnostics(model.06)

##################
# edges + isolates + concurrent
##################
ergm(GeneticNetwork~ edges + isolates + concurrent())
# degenerate

#############################
# edges + concurrent + cycle(3)
#############################
# docs say for undirected graphs cycle value can range from 3 to n
ergm(GeneticNetwork ~ edges + concurrent + cycle(3))
# in addition to normal degeneracy message:
# In addition: Warning messages:
# 1: In mple.existence(pl) : The MPLE does not exist!
#   2: In ergm.mple(nw, fd, m, MPLEtype = MPLEtype, init = init, control = control,  :
#                     glm.fit: fitted probabilities numerically 0 or 1 occurred


#############################
# edges + concurrent() + degree(1) + degree(2) + cycle(3)
#############################
# docs say for undirected graphs cycle value can range from 3 to n
ergm(GeneticNetwork ~ edges + concurrent() + degree(1) + degree(2) + cycle(3))
# degenerate

#################################
# edges + degcor()
#################################
ergm(GeneticNetwork~ edges + degcor())
# degenerate

#################################
# edges + concurrent() + degcor()
#################################
ergm(GeneticNetwork~ edges + concurrent() + degcor())
#degenerate

#######################
# edges + isolatededges
#######################
model.08 = ergm(GeneticNetwork~ edges + isolatededges)
summary(model.08)
plot(gof(model.08))
ComponentGOF(ergm=model.08, 100)
mcmc.diagnostics(model.08)



######################################
# edges + isolatededges + concurrent()
######################################
model.09 = ergm(GeneticNetwork~ edges + isolatededges + concurrent())
summary(model.09)
plot(gof(model.09))
ComponentGOF(ergm=model.09, 100)
mcmc.diagnostics(model.09)



##################
# edges + degree(1) + degree(2) + degree(3) + isolatededges + concurrent()
##################
# summary(GeneticNetwork~edges + degree(1) + degree(2) + degree(3) + isolatededges + concurrent())
ergm(GeneticNetwork~edges + degree(1) + degree(2) + degree(3) + isolatededges + concurrent()) 
# Estimating equations are not within tolerance region.
# Error in ergm.MCMLE(init, nw, model, initialfit = (initialfit <- NULL),  : 
#                       MCMLE estimation stuck. There may be excessive correlation between model terms, suggesting a poor model for the observed data. If target.stats are specified, try increasing SAN parameters.
#                     # isolates? gwesp? gwdegree?


##################
# edges + degree(2) + degree(3) + isolatededges + concurrent()
##################
model.10 = ergm(GeneticNetwork~edges + degree(2) + degree(3) + isolatededges + concurrent())
summary(model.10)
# Call:
#   ergm(formula = GeneticNetwork ~ edges + degree(2) + degree(3) + 
#          isolatededges + concurrent())
# 
# Monte Carlo Maximum Likelihood Results:
#   
#   Estimate Std. Error MCMC % z value Pr(>|z|)    
# edges         -5.28914    0.09061      0 -58.374   <1e-04 ***
#   degree2        0.04286    0.19456      0   0.220   0.8256    
# degree3       -0.34964    0.17166      0  -2.037   0.0417 *  
#   isolatededges -0.19518    0.13104      0  -1.489   0.1364    
# concurrent    -2.49767    0.23848      0 -10.473   <1e-04 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Null Deviance: 1096078  on 790653  degrees of freedom
# Residual Deviance:   15010  on 790648  degrees of freedom
# 
# AIC: 15020  BIC: 15078  (Smaller is better. MC Std. Err. = 1.168)
plot(gof(model.10)) #done
ComponentGOF(ergm=model.10, 100) #done
mcmc.diagnostics(model.10)

##################
# edges + degree(2) + degree(3) + isolatededges + concurrent() + gwesp
##################
#to do model.11 = ergm(GeneticNetwork~ gwesp + edges + degree(2) + degree(3) + isolatededges + concurrent())
summary(model.11) #to do
