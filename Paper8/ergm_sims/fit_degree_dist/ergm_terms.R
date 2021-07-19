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
#15 iterations first time, then 17 the second... what about seed?
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
