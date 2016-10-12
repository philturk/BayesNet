#########SAMPLING FUNCTIONS###########

sample_degree_dist_info <- function(G, num_samples) {

	sample_degree = sample(degree(G), num_samples, replace = TRUE)
	mean_degree_dist_per = tabulate(sample_degree+1)/num_samples
	
	var_degree_dist_per = matrix(0, nrow = length(mean_degree_dist_per), ncol = length(mean_degree_dist_per))
	for (i in c(1:length(mean_degree_dist_per))) {
		for (j in c(1:length(mean_degree_dist_per))) {
			if (i == j) {
				var_degree_dist_per[i,j] = mean_degree_dist_per[i] * (1-mean_degree_dist_per[i]) / num_samples
			} else {
				var_degree_dist_per[i,j] = -mean_degree_dist_per[i] * mean_degree_dist_per[j] / num_samples
			}
		}
	}
	return(list(mean_degree_dist_per, var_degree_dist_per))
}

var_deg_dist_compute <- function(g_deg_dist_orig) {

	sample_size = sum(g_deg_dist_orig)

	p = g_deg_dist_orig / sum(g_deg_dist_orig)

	var_x = matrix(data = 0, nrow = length(p), ncol = length(p))

	for (i in c(1:length(p))) {
		for (j in c(1:i)) {
			if (i == j) {
				var_x[i,j] = p[i]*(1-p[i])
			} else {
				var_x[j,i] = var_x[i,j] = -p[i]*p[j]
			}
		}
	}

#	var_x = 1/sample_size * var_x

	return(var_x)

}

#Check if necessary

sample_degdist_g <- function(g = NULL, mean_deg_dist, num_node_samples, deg_corr = TRUE) {

	if (deg_corr == TRUE) {
		mean_deg_dist_adj = mean_deg_dist*num_node_samples
		zero_deg_dist_loc = which(mean_deg_dist_adj == 0)	
	
		mean_deg_dist_adj[zero_deg_dist_loc] = mean_deg_dist_adj[zero_deg_dist_loc] + 1
		mean_deg_dist = mean_deg_dist_adj / sum(mean_deg_dist_adj)
	}

	var_x = matrix(0, nrow = length(mean_deg_dist), ncol = length(mean_deg_dist))
	for (i in c(1:(length(mean_deg_dist)))) {
		for (j in c(1:(length(mean_deg_dist)))) {
			if (i == j) {
				var_x[i,j] = mean_deg_dist[i]*(1-mean_deg_dist[i]) / num_node_samples
			} else {
				var_x[i,j] = -mean_deg_dist[i]*mean_deg_dist[j] / num_node_samples
			}
		}
	}
	return(list(var_x, mean_deg_dist))
}

sample_mixing_g <- function(g = NULL, mean_mixing, num_mixing_samples) {

	mean_mixing = mean_mixing[upper.tri(mean_mixing, diag = TRUE)]

	var_x = matrix(0, nrow = length(mean_mixing), ncol = length(mean_mixing))
	for (i in c(1:(length(mean_mixing)))) {
		for (j in c(1:(length(mean_mixing)))) {
			if (i == j) {
				var_x[i,j] = mean_mixing[i]*(1-mean_mixing[i]) / num_mixing_samples
			} else {
				var_x[i,j] = -mean_mixing[i]*mean_mixing[j] / num_mixing_samples
			}
		}
	}
	return(var_x)
}



