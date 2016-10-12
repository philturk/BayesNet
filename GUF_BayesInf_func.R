
BI_posterior_prior <- function(g, Network_stats, Prob_Distr_Params) {

	g_orig = g

	if  (Network_stats == 'Density') {

	  mean_density = Prob_Distr_Params[[1]][[1]]
	  var_density = Prob_Distr_Params[[1]][[2]]
	  
		g_density_orig =  network.edgecount(g_orig) / choose(network.size(g_orig),2)
		var_density_orig = (g_density_orig * (1-g_density_orig))

		mean_density_1 = ((mean_density/var_density) + (network.edgecount(g_orig)/var_density_orig))
		mean_density_2 = ((1/var_density) + (choose(network.size(g_orig),2)/var_density_orig))
		mean_density = mean_density_1 / mean_density_2

		var_density = 1/mean_density_2

		return(list(list(mean_density, var_density)))
	}

	if (Network_stats == 'DegreeDist') {

	  mean_deg_dist = Prob_Distr_Params[[1]][[1]]/population
	  var_deg_dist = Prob_Distr_Params[[1]][[2]]/population^2

	#Likelihood is normal should be multinominal (close so not a big deal)

		mean_deg_dist_orig = tabulate(degree(g_orig)/2 + 1) 
		mean_deg_dist_orig = c(mean_deg_dist_orig,rep(0,max(0,length(mean_deg_dist)-length(mean_deg_dist_orig))))
		mean_deg_dist_orig = (mean_deg_dist_orig + .2) / sum(mean_deg_dist_orig + .2) #####THE +0.2 is to prevent singular matrices for variance

		var_deg_dist_orig = var_deg_dist_compute(mean_deg_dist_orig)
		inv_var_deg_dist_orig = solve(var_deg_dist_orig[-length(mean_deg_dist_orig), -length(mean_deg_dist_orig)]) #NEW

		inv_var_deg_dist = solve(var_deg_dist[-length(mean_deg_dist), -length(mean_deg_dist)]) #NEW

		mean_deg_dist_1 = ((inv_var_deg_dist %*% (mean_deg_dist[-length(mean_deg_dist)])) + (network.size(g) * (inv_var_deg_dist_orig %*% (mean_deg_dist_orig[-length(mean_deg_dist)]))))
		mean_deg_dist_2 = solve(inv_var_deg_dist + (network.size(g) * inv_var_deg_dist_orig)) #NEW
		mean_deg_dist = mean_deg_dist_2 %*% mean_deg_dist_1
		var_deg_dist = mean_deg_dist_2

	#Correct for removed last row / column

		mean_deg_dist = c(mean_deg_dist, 1-sum(mean_deg_dist))
		sample_size_var = (mean_deg_dist[1] * (1-mean_deg_dist[1])) / var_deg_dist[1,1]	

		var_deg_dist_lastrow = array(dim = length(mean_deg_dist)-1)
		for (var_counter in c(1:(length(mean_deg_dist)-1))) {
			var_deg_dist_lastrow[var_counter] = -mean_deg_dist[length(mean_deg_dist)] * mean_deg_dist[var_counter] / sample_size_var
		}
		var_deg_dist = rbind(var_deg_dist, var_deg_dist_lastrow)
		p_last = mean_deg_dist[length(mean_deg_dist)]
		var_deg_dist = cbind(var_deg_dist, c(var_deg_dist_lastrow, (p_last * (1-p_last))/ sample_size_var))

#print(c("MEAN: ", mean_deg_dist))
#print(c("VAR : ", diag(var_deg_dist)))
		mean_deg_dist = mean_deg_dist * population
		var_deg_dist = var_deg_dist * (population^2)
		  
		return(list(list(mean_deg_dist, var_deg_dist)))

	}

	if (Network_stats == 'DegMixing') {

	  mean_dmm = Prob_Distr_Params[[1]][[1]]
	  var_dmm = Prob_Distr_Params[[1]][[2]]

	#Likelihood is normal should be multinominal (close so not a big deal)

		mean_dmm_orig = g_dmm

		mean_dmm_orig = cbind(mean_dmm_orig, matrix(0, nrow = dim(mean_dmm_orig)[1], ncol = dim(mean_dmm)[2] - dim(mean_dmm_orig)[2]))
		mean_dmm_orig = rbind(mean_dmm_orig, matrix(0, nrow = dim(mean_dmm)[1] - dim(mean_dmm_orig)[1], ncol = dim(mean_dmm_orig)[2]))
		mean_dmm_orig = (mean_dmm_orig + .2)/ sum(mean_dmm_orig[upper.tri(mean_dmm_orig, diag = TRUE)] + .2) #####THE +0.2 is to prevent singular matrices for variance

		mean_dmm_orig_vector = mean_dmm_orig[upper.tri(mean_dmm_orig, diag = TRUE)]
		var_dmm_orig = var_deg_dist_compute(mean_dmm_orig_vector)
		inv_var_dmm_orig = solve(var_dmm_orig[-length(mean_dmm_orig_vector), -length(mean_dmm_orig_vector)]) #NEW

		mean_dmm_vector = mean_dmm[upper.tri(mean_dmm, diag = TRUE)]
		inv_var_dmm = solve(var_dmm[-length(mean_dmm_vector), -length(mean_dmm_vector)])

		mean_dmm_1 = ((inv_var_dmm %*% (mean_dmm_vector[-length(mean_dmm_vector)])) + (network.edgecount(g) * (inv_var_dmm_orig %*% (mean_dmm_orig_vector[-length(mean_dmm_orig_vector)]))))
		mean_dmm_2 = solve(inv_var_dmm + (network.edgecount(g) * inv_var_dmm_orig)) #NEW
		mean_dmm_vector = mean_dmm_2 %*% mean_dmm_1
		var_dmm = mean_dmm_2

	#Correct for removed last row / column

		mean_dmm_vector = c(mean_dmm_vector, 1-sum(mean_dmm_vector))
		sample_size_var = (mean_dmm_vector[1] * (1-mean_dmm_vector[1])) / var_dmm[1,1]	

		var_dmm_lastrow = array(dim = length(mean_dmm_vector)-1)
		for (var_counter in c(1:(length(mean_dmm_vector)-1))) {
			var_dmm_lastrow[var_counter] = -mean_dmm_vector[length(mean_dmm_vector)] * mean_dmm_vector[var_counter] / sample_size_var
		}
		var_dmm = rbind(var_dmm, var_dmm_lastrow)
		p_last = mean_dmm_vector[length(mean_dmm_vector)]
		var_dmm = cbind(var_dmm, c(var_dmm_lastrow, (p_last * (1-p_last))/ sample_size_var))

		mean_dmm = matrix(nrow = dim(mean_dmm)[1], ncol = dim(mean_dmm)[2])
		mean_dmm[upper.tri(mean_dmm, diag = TRUE)] = mean_dmm_vector
		mean_dmm[lower.tri(mean_dmm, diag = TRUE)] = t(mean_dmm)[lower.tri(mean_dmm, diag = TRUE)]

#print(c("MEAN: ", mean_dmm))
#print(c("VAR : ", diag(var_dmm)))

		return(list(list(mean_dmm, var_dmm)))

	}

}


#for the probability of accept of reject

update_r <- function(r, g, P, Ia, Il, R_times, edge1, beta_a, beta_l) {

	add_edge = 1
	if (is.adjacent(g, edge1[1], edge1[2])) { 		#Remove edge
		add_edge = 0
	} 

	if (add_edge == 1) {
		p_edge1 = r / (1 + r)
	} else {
		p_edge1 = 1 / (1 + r)
	}


	if ((is.adjacent(P, edge1[1], edge1[2])) || (is.adjacent(P, edge1[2], edge1[1])) ) {
		return(0)
	} else {
		if ((Ia[edge1[1]] == Inf) && (Ia[edge1[2]] == Inf)) {
			return(r)
		} else {
			Il_i = Il[edge1[1]]
			Ia_i = Ia[edge1[1]]
			R_i = R[edge1[1]]
			Il_j = Il[edge1[2]]
			Ia_j = Ia[edge1[2]]
			R_j = R[edge1[2]]

			if (Ia[edge1[2]] < Ia[edge1[1]]) {
				time_a = min(Il_j,Ia_i)-Ia_j
				time_l = max(min(R_j,Ia_i),Il_j) - Il_j
				muij = (1-pexp(time_a,rate = beta_a))*(1-pexp(time_l, rate = beta_l))
			} else {
				time_a = min(Il_i,Ia_j)-Ia_i
				time_l = max(min(R_i,Ia_j),Il_i) - Il_i
				muij = (1-pexp(time_a,rate = beta_a))*(1-pexp(time_l, rate = beta_l))
			}
			p_noinfect = (muij*p_edge1)/((1-p_edge1) + muij*p_edge1)

			if (add_edge == 1) {
				return(p_noinfect / (1 - p_noinfect) )
			} else {
				return((1-p_noinfect) / p_noinfect )
			}
		}
	}
}



