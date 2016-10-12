G.generate <- function(g_type = 2, population = 100, bool_ER = TRUE, bool_Assort = FALSE, bool_Clustering = FALSE, ER_prob = .05, Assort_val = 0, Cluster_val = 0) {

	if (g_type == 1) {
		G = network.initialize(population, directed = FALSE)
		G = add.edges(G, c(0,0,0,1,3,3,4,4,4,8)+1, c(1,3,5,2,4,8,5,6,7,9)+1)
	}
	if (g_type == 2) {
		if (bool_Clustering == TRUE) {

			library('ergm')
			edges = runif(population*(population-1)/2) < ER_prob
			g.use <- network.initialize(population,directed=FALSE)
			g.use[upper.tri(g.use)] = edges

			g.sim <- simulate(~gwesp(0.2,fixed=TRUE), nsim=1, coef=c(Cluster_val),basis=g.use, constraints = ~ degreedist)

			y_edges_c = as.matrix.network(g.sim,matrix.type="edgelist", as.sna.edgelist = TRUE)
			y_edges_c = y_edges_c[c(1:(dim(y_edges_c)[1]/2)),c(1:2)	]
				
			G = network.initialize(population, directed = FALSE)
			G = add.edges(G, y_edges_c[,1], y_edges_c[,2])

		}
		if (bool_ER == TRUE) {
			edges = runif(population*(population-1)/2) < ER_prob
			G <- network.initialize(population,directed=FALSE)
			G[upper.tri(G)] = edges
		}
		if (bool_Assort == TRUE) {
			library('ergm')
			edges = runif(population*(population-1)/2) < ER_prob
			G <- network.initialize(population,directed=FALSE)
			G[upper.tri(G)] = edges
			if (abs(Assort_val) > 0) {
				G = ergm(G~degcor, constraints = ~degreedist, target.stats = Assort_val)$newnetwork
			}
		}
	}
	return(G)
}



load_data <- function(cluster_option) {
	if (cluster_option == 1) {
		setwd("~/Data")
	} else {
		setwd("C:\\Users\\Ravi\\Desktop\\Ravi\\Network Research\\Helleringer\\LNS-Data") # 
	}
	edges <- read.table("LNS_arc.csv", header = TRUE, sep = ",")

	#edges = edges[(which(edges[,"tag_edge"] == 0)),] # only take edges that both declared
	#edges = edges[(which(edges[,2] == "Female")),]

	#All the to_edges are male and all the from edges are female
	#There are some issues with some nodes
	#intersect(to_edges, from_edges) #=  781016 100102 120121 140030
	#edges[which(edges[,1] == 781016),]
	#edges[which(edges[,1] == 100102),]
	#edges[which(edges[,1] == 120121),]
	#edges[which(edges[,3] == 120121),]
	#edges[which(edges[,1] == 140030),]
	#edges[which(edges[,3] == 140030),]

	error_edges = union(which(edges[,1] == 781016),which(edges[,1] == 100102))
	error_edges = union(error_edges,which(edges[,3] == 100102))
	error_edges = union(error_edges,which(edges[,1] == 120121))
	error_edges = union(error_edges,which(edges[,3] == 120121))
	error_edges = union(error_edges,which(edges[,1] == 140030))
	error_edges = union(error_edges,which(edges[,3] == 140030))
	error_edges = unique(error_edges)

	edges = edges[-error_edges,]
	edges
}

generate_sexual_network <- function(edges, edge_types) {

	if (edge_types == 1) {
		edges1 <- which(edges[,2] == "Female")
	}	
	if (edge_types == 2) {
		edges1 <- which(edges[,2] == "Male")
	}	
	if (edge_types == 3) {
		num_edges = length(edges[,1])
		duplicate_edges = c()
		original_edges = c()
		for (i in c(1:(num_edges-1))) {
			set1 = which(edges[c((i+1):num_edges),1] == edges[i,1]) + i
			set2 = which(edges[c((i+1):num_edges),3] == edges[i,1]) + i
 			set3 = union(set1,set2)
			for (j in set3) {
				if ((edges[i,1] == edges[j,1]) && (edges[i,3] == edges[j,3])) {
					duplicate_edges = c(duplicate_edges, j)
					original_edges = c(original_edges, i)
				}
				if ((edges[i,3] == edges[j,1]) && (edges[i,1] == edges[j,3])) {
					duplicate_edges = c(duplicate_edges, j)
					original_edges = c(original_edges, i)
				}
			}
		}
		edges <- c(original_edges,duplicate_edges)
	}


	from_edges <- edges[,3]
	gender <- edges[,2]
	to_edges <- edges[,1]

	for (i in 1:length(gender)) {
		if (gender[i] == "Female") {
			temp = from_edges[i] 
			from_edges[i] = to_edges[i]
			to_edges[i] = temp
		}	
	}

	#intersect(to_edges, from_edges) #= NULL which what we want

	network_edges <- c()
	for (i in 1:length(from_edges)) {
		network_edges <- c(network_edges, from_edges[i], to_edges[i])
	}

	#Network_edges has female than male listed....from_edges are female

	u_network_edges <- unique(network_edges)
	network_edge_orig <- network_edges

	network_edges2 <- network_edges

	for (i in 1:length(network_edges)){
		new_id <- which(network_edges[i] == u_network_edges)
		network_edges[i] <- new_id - 1	
	}

	pop <- max(network_edges+1)

	g <- graph.empty(n = pop, directed=FALSE)

	if (edge_types > 0) {
		network_edges1 = network_edges
		dim(network_edges1) = c(2, length(network_edges1)/2)
		network_edges1 = t(network_edges1)
		network_edges1 = network_edges1[edges1,]
		network_edges1 = c(network_edges1)
		g <- add.edges(g, network_edges1)
	} else {
		g <- add.edges(g, network_edges)
	}
	g <- simplify(g)

	dim(network_edges) = c(2,length(network_edges)/2)
	females = unique(network_edges[1,])
	males = unique(network_edges[2,])

	g <- set.vertex.attribute(g, "Gender", index=females, "Female")
	g <- set.vertex.attribute(g, "Gender", index=males, "Male")

	for (i in c(0:(pop-1))) {
		g <- set.vertex.attribute(g, "Degree", index=i, value = length(neighborhood(g, 1, i)[[1]])-1)
	}

	g
}


generate_graph_ex <- function(population, gen_g_option, deg_zero_num, max_d) {

	if (gen_g_option == 0) {
		deg_g = c(10,30,16)

		g = degree.sequence.game(c(rep(1,deg_g[2]),rep(2,deg_g[3])))
		g = simplify(g)
		g = add.vertices(g, deg_g[1])  #Number of degree 0 nodes
		g = set.vertex.attribute(g, "race", index=V(g), value = 1)
		g = set.vertex.attribute(g, "race", index=(sample(c(1:length(V(g))),length(V(g))/2,replace=FALSE)-1), value = 2)
	}
	if (gen_g_option == 1) {
		graph_file = paste("Collection_graphs_", population, "_022812.RData", sep="") #File Name of generated networks
		setwd("C:\\Users\\Ravi\\Dropbox\\Ravi Network\\Concurrency")
		load(graph_file)
	
		graph_id = sample(c(1:100),1)
		g = collection_g[[graph_id]]

	print(graph_id)

		g = add.vertices(g,deg_zero_num)

		g = set.vertex.attribute(g, "race", index=V(g), value = 1)
		g = set.vertex.attribute(g, "race", index=(sample(c(1:length(V(g))),length(V(g))/2,replace=FALSE)-1), value = 2)
	} 
	if (gen_g_option == 2) {
		g = generate_graph(population, degree_type="sim", max_d = max_d)
		g = add.vertices(g,deg_zero_num)
		g = set.vertex.attribute(g, "race", index=V(g), value = 1)
		g = set.vertex.attribute(g, "race", index=(sample(c(1:length(V(g))),length(V(g))/2,replace=FALSE)-1), value = 2)
	}
	g = simplify(g)
	return(g)
}


add_health_data <- function(population, data_dir, prob_type, multi_school = NULL, school_loc = NULL) {

  
  network.data = read.dta(paste(data_dir, "21600-0018-Data.dta", sep = ""))
  individual.data = read.dta(paste(data_dir, "21600-0001-Data.dta", sep = ""))
  set2.data = read.dta(paste(data_dir, "21600-0002-Data.dta", sep = ""))
  
	core_sample = individual.data[,"SMP01"]

	network.data = network.data[which(core_sample == "Yes"),]
	individual.data = individual.data[which(core_sample == "Yes"),]
	set2.data = set2.data[which(core_sample == "Yes"),]

	gender = individual.data[,"BIO_SEX"]
	school_id = unique(set2.data[,2])
	Num_female = round(network.data[,"NASS2"] * network.data[,"AXSS2"])
	Num_male = network.data[,"NASS2"] - Num_female
	node_ID = network.data[,"AID"]

	dim(network.data)

	gender_mixing = matrix(0, nrow = 3, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
	gender_mixing_norm = matrix(0, nrow = 3, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF

	degree_dist = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
	degree_dist_norm = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF

	degree_dist_F = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
	degree_dist_norm_F = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF

	degree_dist_M = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
	degree_dist_norm_M = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF

	counter = 1

	school_lengths = c()

	for (i in school_id) {
		same_school = which(set2.data[,2] == i)
	
		school_lengths = c(school_lengths, length(same_school))

		for (j in same_school) {
			gender_j = gender[j]
			num_female_j = if(is.na(Num_female[j])){0}else{Num_female[j]}
			num_male_j = if(is.na(Num_male[j])){0}else{Num_male[j]}

			if (gender_j == "Male") {
				gender_mixing[1,counter] = gender_mixing[1,counter] + num_male_j
				gender_mixing[2,counter] = gender_mixing[2,counter] + num_female_j
			}
			if (gender_j == "Female") {
				gender_mixing[3,counter] = gender_mixing[3,counter] + num_female_j
				gender_mixing[2,counter] = gender_mixing[2,counter] + num_male_j	
			}
			if (!is.na(network.data[j,"NASS2"])) {
				degree_dist[(network.data[j,"NASS2"]+1),counter] = degree_dist[(network.data[j,"NASS2"]+1),counter] + 1
				if (gender_j == "Female") {
					degree_dist_F[(network.data[j,"NASS2"]+1),counter] = degree_dist_F[(network.data[j,"NASS2"]+1),counter] + 1
				}
				if (gender_j == "Male") {
					degree_dist_M[(network.data[j,"NASS2"]+1),counter] = degree_dist_M[(network.data[j,"NASS2"]+1),counter] + 1
				}
			}
		}
		degree_dist_norm[,counter] = 	degree_dist[,counter]/sum(degree_dist[,counter])
		degree_dist_norm_F[,counter] = degree_dist_F[,counter]/sum(degree_dist_F[,counter])
		degree_dist_norm_M[,counter] = degree_dist_M[,counter]/sum(degree_dist_M[,counter])

		gender_mixing_norm[,counter] = gender_mixing[,counter]/sum(gender_mixing[,counter])

		counter = counter + 1
	}

	if (multi_school == FALSE) {

		num_males = sum(degree_dist_M[,school_loc])
		num_females = sum(degree_dist_F[,school_loc])
		num_mixing = sum(gender_mixing[,school_loc])

		mean_mixing = gender_mixing_norm[,school_loc]
		mean_deg_dist = degree_dist_norm[,school_loc]
		mean_deg_dist_F = degree_dist_norm_F[,school_loc]
		mean_deg_dist_M = degree_dist_norm_M[,school_loc]

		while(mean_deg_dist[length(mean_deg_dist)] == 0) {
			mean_deg_dist = mean_deg_dist[-length(mean_deg_dist)]
		}
		mean_deg_dist_F = mean_deg_dist_F[c(1:length(mean_deg_dist))]
		mean_deg_dist_M = mean_deg_dist_M[c(1:length(mean_deg_dist))]

		mean_mixing_matrix = matrix(0, nrow = 2, ncol = 2)
		mean_mixing_matrix[1,1] = mean_mixing[1]
		mean_mixing_matrix[2,1] = mean_mixing[2]
		mean_mixing_matrix[1,2] = mean_mixing[2]
		mean_mixing_matrix[2,2] = mean_mixing[3]
	
	} else {
		mean_mixing = apply(gender_mixing_norm, 1, mean, na.rm=TRUE)
		mean_deg_dist = apply(degree_dist_norm, 1, mean, na.rm=TRUE)

		mean_deg_dist_F = NULL
		mean_deg_dist_M = NULL

		num_males = NULL
		num_females = NULL
		num_mixing = NULL
	}

	degree_dist_g = c()
	for (i in c(1:length(mean_deg_dist))) {
		degree_dist_g = c(degree_dist_g, rep((i-1), trunc(population*mean_deg_dist[i])))
	}

	if (multi_school == FALSE) {
		if (school_loc == 78) {
			if (population-length(degree_dist_g) > 0) {
				degree_dist_g = c(degree_dist_g, c(1,1,1,2))
			}
		}
		if (sum(school_loc == c(4,84,90)) > 0) {
			if (population-length(degree_dist_g) > 0) {
				degree_dist_g = c(degree_dist_g, rep(5,population-length(degree_dist_g)))
			}
		}
	} else {
		if (population-length(degree_dist_g) > 0) {
			degree_dist_g = c(degree_dist_g, rep(5,population-length(degree_dist_g)))
		}
	}

	g = degree.sequence.game(degree_dist_g)
	g = simplify(g)

	g = set.vertex.attribute(g, "race", index=V(g), value=1)
	if (multi_school == TRUE) {
		g = set.vertex.attribute(g, "race", index=sample(V(g),(population/2), replace=FALSE), value=2)
	} else {
		num_females_a = round(population * (sum(degree_dist_F[,school_loc])/sum(degree_dist[,school_loc])))
		g = set.vertex.attribute(g, "race", index=(sample((V(g)+1),num_females_a, replace=FALSE)-1), value=2)
	}
			
	mean_mixing_matrix = matrix(0, nrow = 2, ncol = 2)
	mean_mixing_matrix[1,1] = mean_mixing[1]
	mean_mixing_matrix[2,1] = mean_mixing[2]
	mean_mixing_matrix[1,2] = mean_mixing[2]
	mean_mixing_matrix[2,2] = mean_mixing[3]

	if (prob_type == 5) {  #Kernel Density	
		remove_id = which(is.nan(degree_dist_norm[1,]))
		bw <- npudensbw(dat=t(degree_dist_norm[,-remove_id]))
	} else {  #Multivariate Density
		bw = NULL
	}

	if (prob_type == 6) {  #Dirichlet Density	
		remove_id = which(is.nan(degree_dist_norm[1,]))

		degree_dist_norm2 = degree_dist_norm[,-remove_id]
		degree_dist_norm2 = c(degree_dist_norm2)
		degree_dist_norm2[which(degree_dist_norm2 == 0)] = runif(length(which(degree_dist_norm2 == 0)), min=.00001,max = .0001)
				
		dim(degree_dist_norm2) = c(11,length(degree_dist_norm2)/11)

		for (i in c(1:dim(degree_dist_norm2)[2])) {
			degree_dist_norm2[,i] = degree_dist_norm2[,i] / sum(degree_dist_norm2[,i])
		}
		
		valid_length = which(school_lengths >= 50)
		
		dirichlet_dist <- exp(coef(vglm(t(degree_dist_norm2)~ 1, dirichlet, trace = FALSE, crit="c"))) #[-1,valid_length]
	} else {  #Multivariate Density
		dirichlet_dist = NULL
	}

	return(list(g, mean_deg_dist,mean_mixing_matrix, bw, dirichlet_dist, mean_deg_dist_M, mean_deg_dist_F, num_males, num_females, num_mixing))

}



add_health_data_2 <- function(data_dir=NULL) {
  
  network.data = read.dta(paste(data_dir, "21600-0018-Data.dta", sep = ""))
  individual.data = read.dta(paste(data_dir, "21600-0001-Data.dta", sep = ""))
  set2.data = read.dta(paste(data_dir, "21600-0002-Data.dta", sep = ""))
  
  core_sample = individual.data[,"SMP01"]
  
  network.data = network.data[which(core_sample == "Yes"),]
  individual.data = individual.data[which(core_sample == "Yes"),]
  set2.data = set2.data[which(core_sample == "Yes"),]
  
  gender = individual.data[,"BIO_SEX"]
  school_id = unique(set2.data[,2])
#  Num_female = round(network.data[,"NASS2"] * network.data[,"AXSS2"])
#  Num_male = network.data[,"NASS2"] - Num_female

  Num_female = round(network.data[,"NAS2"] * network.data[,"AXS2"])
  Num_male = network.data[,"NAS2"] - Num_female
  
  node_ID = network.data[,"AID"]
  
  gender_mixing = matrix(0, nrow = 3, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
  gender_mixing_norm = matrix(0, nrow = 3, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
  
#  degree_dist = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
#  degree_dist_norm = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
  
#  degree_dist_F = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
#  degree_dist_norm_F = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
  
#  degree_dist_M = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
#  degree_dist_norm_M = matrix(0, nrow = 11, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF

  degree_dist = matrix(0, nrow = 33, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
  degree_dist_norm = matrix(0, nrow = 33, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
  
  degree_dist_F = matrix(0, nrow = 33, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
  degree_dist_norm_F = matrix(0, nrow = 33, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
  
  degree_dist_M = matrix(0, nrow = 33, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
  degree_dist_norm_M = matrix(0, nrow = 33, ncol = length(school_id)) #row 1 - MM, row 2 - MF, row 3 - FF
  
  counter = 1
  
  school_lengths = c()
  
  for (i in school_id) {
    same_school = which(set2.data[,2] == i)
    
    school_lengths = c(school_lengths, length(same_school))
    
    for (j in same_school) {
      gender_j = gender[j]
      num_female_j = if(is.na(Num_female[j])){0}else{Num_female[j]}
      num_male_j = if(is.na(Num_male[j])){0}else{Num_male[j]}
      
      if (gender_j == "Male") {
        gender_mixing[1,counter] = gender_mixing[1,counter] + num_male_j
        gender_mixing[2,counter] = gender_mixing[2,counter] + num_female_j
      }
      if (gender_j == "Female") {
        gender_mixing[3,counter] = gender_mixing[3,counter] + num_female_j
        gender_mixing[2,counter] = gender_mixing[2,counter] + num_male_j  
      }
#      if (!is.na(network.data[j,"NASS2"])) {
#        degree_dist[(network.data[j,"NASS2"]+1),counter] = degree_dist[(network.data[j,"NASS2"]+1),counter] + 1
#        if (gender_j == "Female") {
#          degree_dist_F[(network.data[j,"NASS2"]+1),counter] = degree_dist_F[(network.data[j,"NASS2"]+1),counter] + 1
#        }
#        if (gender_j == "Male") {
#          degree_dist_M[(network.data[j,"NASS2"]+1),counter] = degree_dist_M[(network.data[j,"NASS2"]+1),counter] + 1
#        }
#      }
      
      if (!is.na(network.data[j,"NAS2"])) {
        degree_dist[(network.data[j,"NAS2"]+1),counter] = degree_dist[(network.data[j,"NAS2"]+1),counter] + 1
        if (gender_j == "Female") {
          degree_dist_F[(network.data[j,"NAS2"]+1),counter] = degree_dist_F[(network.data[j,"NAS2"]+1),counter] + 1
        }
        if (gender_j == "Male") {
          degree_dist_M[(network.data[j,"NAS2"]+1),counter] = degree_dist_M[(network.data[j,"NAS2"]+1),counter] + 1
        }
      }
      
    }
    degree_dist_norm[,counter] = 	degree_dist[,counter]/sum(degree_dist[,counter])
    degree_dist_norm_F[,counter] = degree_dist_F[,counter]/sum(degree_dist_F[,counter])
    degree_dist_norm_M[,counter] = degree_dist_M[,counter]/sum(degree_dist_M[,counter])
    
    gender_mixing_norm[,counter] = gender_mixing[,counter]/sum(gender_mixing[,counter])
    
    counter = counter + 1
  }
  
  remove_id = unique(c(which(is.nan(degree_dist_norm[1,])), which(is.nan(degree_dist_norm_F[1,])), which(is.nan(degree_dist_norm_M[1,])), which(is.nan( gender_mixing_norm[1,]))))
  
  degree_dist_norm = degree_dist_norm[,-remove_id]
  degree_dist_norm_F = degree_dist_norm_F[,-remove_id]
  degree_dist_norm_M = degree_dist_norm_M[,-remove_id]
  gender_mixing_norm = gender_mixing_norm[,-remove_id]
  
  return(list(degree_dist_norm, degree_dist_norm_F, degree_dist_norm_M, gender_mixing_norm))
  
}



initial_bipartite_community_graph <- function(num_communities, populations) {

	pop1_f = populations[1] / 2
	pop1_m = populations[1] / 2
	pop2_f = populations[2] / 2
	pop2_m = populations[2] / 2

	pop1 = pop1_f + pop1_m
	pop2 = pop2_f + pop2_m

	g <- graph.empty(n = (pop1_f + pop1_m + pop2_f + pop2_m), directed=FALSE)

	row.sums <- simnb(n=pop1_f, v=c(5, .7), maxdeg=min(pop1_m-1, 7))
	col.sums <- sample(row.sums, pop1_m)

	bg <- network.form(row.sums, col.sums)
	bge = c(unlist(bg$mel))
	dim(bge) = c(3, length(bge)/3)
	bge = bge[c(1,2),]
	bge = c(bge)
	gender_a1 = c(rep("Female",length(row.sums)), rep("Male",length(col.sums)))

	g <- add.edges(g, (bge-1))

	#row.sums <- simnb(n=pop2_f, v=c(5, .7), maxdeg=max(pop2_m, 7))
	#col.sums <- sample(row.sums, pop2_m)
	bg <- network.form(row.sums, col.sums)
	bge = c(unlist(bg$mel))
	dim(bge) = c(3, length(bge)/3)
	bge = bge[c(1,2),]
	bge = c(bge)
	gender_a2 = c(rep("Female",length(row.sums)), rep("Male",length(col.sums)))

	g <- add.edges(g, (bge-1+pop1))

	gender_a = c(gender_a1, gender_a2)	
	communities = c(rep(0,pop1), rep(1,pop2))

	g = set.vertex.attribute(g, "Gender", value = gender_a)
	g <- set.vertex.attribute(g, "Degree", value = degree(g))
	g <- set.vertex.attribute(g, "Community", value = communities)

	g <- simplify(g)

	return(g)
}



