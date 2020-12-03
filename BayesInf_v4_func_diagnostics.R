#######Diagnostics###########

genetic_seq_diagnostic <- function(P,T) {
  
  dist_linked_cases = c()
  for (i in c(1:ecount(P))) {
    dist_linked_cases = c(dist_linked_cases, T[(P[[3]][i])+1, (P[[4]][i])+1])
  }
  
  d1 = density(dist_linked_cases)
  
  T_data = c(T[upper.tri(T,diag = FALSE)])
  d2 = density(T_data[which(T_data > 0)])
  
  plot(c(min(d1$x,d2$x),max(d1$x,d2$x)),c(0,max(d1$y,d2$y)), type = "n", main = "Linked vs General", xlab = "Dist in Phlyogenetics", ylab="Density")
  
  polygon(d1, col = '#00FF3348')
  polygon(d2, col = '#00003348')
  
  return(list(mean(dist_linked_cases),mean(T_data[which(T_data > 0)]), 
              var(dist_linked_cases),var(T_data[which(T_data > 0)])))
}

uncategorized_diagnostic <- function(P_a, G_a, G, P, Ia, Il, R, Initial_Data) {
  
  print(cor(Ia, Initial_Data[[3]]))
  print(cor(Il, Initial_Data[[4]]))
  
}

accuracy_G_P_diagnostic <- function(P_a, G_a, Initial_Data) {
  
  G = Initial_Data[[1]]
  P = Initial_Data[[2]]
  
  adj_m = get.adjacency(P)
  for (i in c(51:length(P_a))) {
    adj_m  = adj_m + get.adjacency(P_a[[i]])
  }
  adj_m = adj_m / length(P_a)
  adj_P = get.adjacency(P_a[[1]])
  
  false_positives_dist_P = adj_m[which(adj_P == 0)]
  true_positives_dist_P = adj_m[which(adj_P == 1)]
  
  par(mfrow = c(1,2))
  
  d1 = density(false_positives_dist_P, bw = .01)
  d2 = density(true_positives_dist_P, bw = .01)
  
  plot(c(0,1), c(0,40), type = "n", main = "Transmission Network", xlab = "Percent Edge Occurred", ylab="Density")
  
  polygon(d1, col = '#00FF3348')
  polygon(d2, col = '#00003348')
  
  print(mean(false_positives_dist_P))
  print(mean(true_positives_dist_P))
  
  adj_Gm = get.adjacency(G_a[[2]])
  for (i in c(3:length(G_a))) {
    adj_Gm  = adj_Gm + get.adjacency(G_a[[i]])
  }
  adj_Gm = adj_Gm / length(G_a)
  adj_G = get.adjacency(G_a[[1]])
  
  false_positives_dist_G = adj_Gm[which(adj_G == 0)]
  true_positives_dist_G = adj_Gm[which(adj_G == 1)]
  
  d1 = density(false_positives_dist_G, bw = .01)
  d2 = density(true_positives_dist_G, bw = .01)
  
  plot(c(0,1), c(0,40), type = "n", main = "Contact Network", xlab = "Percent Edge Occurred", ylab="Density")
  
  polygon(d1, col = '#00FF3348')
  polygon(d2, col = '#00003348')
  
  print(mean(false_positives_dist_G))
  print(mean(true_positives_dist_G))
  
}


network_properties_diagnostic <- function(G_a, burnin, view_graphs = FALSE) {
  
  burnin = burnin + 1
  population = length(network.size(G_a[[1]]))
  
  ecount_truth = network.edgecount(G_a[[1]])
  assort_truth =  summary(G_a[[1]]~degcor)
  degree_dist_truth = tabulate(degree(G_a[[1]])/2+1)
  
  ecount_a = c()
  for (i in c(burnin:length(G_a))) {
    ecount_a = c(ecount_a, network.edgecount(G_a[[i]]))
  }
  
  yrange = c( .99*min(c(cumsum(ecount_a) / seq_along(ecount_a), ecount_truth)) , 1.01*max(c(cumsum(ecount_a) / seq_along(ecount_a), ecount_truth)))
  
  if (view_graphs == TRUE) {
    plot(c(1:length(ecount_a)),cumsum(ecount_a) / seq_along(ecount_a), main = "Convergence of Density", ylab = "Cum. Mean", xlab = "Iterations", ylim = yrange )
    abline(h=ecount_truth, col = "red")
    readline(prompt = "Pause. Press <Enter> to continue...")
    
    plot_data = ecount_a/choose(network.size(G_a[[1]]),2)
    plot(stats::density(plot_data), main = "Density plot Network Property Edge Density", xlab = "Edge Density")
    abline(v=(ecount_truth/choose(network.size(G_a[[1]]),2)), col = "red") 
    readline(prompt = "Pause. Press <Enter> to continue...")
  }
  
  max_deg = 0
  for (i in c(burnin:length(G_a))) {
    if (max_deg < max(degree(G_a[[i]])/2)) {
      max_deg = max(degree(G_a[[i]])/2)
    }
  }
  max_deg = max(max_deg, max(degree(G_a[[1]])/2))
  
  assort_a = c()
  degree_dist_mat = matrix(0, nrow = length(c(burnin:length(G_a))), ncol = max_deg+1)
  for (i in c(burnin:length(G_a))) {
    index_i = i - burnin + 1
    assort_a = c(assort_a, summary(G_a[[i]]~degcor))
    degree_dist_mat[index_i,c(1:length(tabulate(degree(G_a[[i]])/2+1)))] = tabulate(degree(G_a[[i]])/2+1)
  }
  
  if (view_graphs == TRUE) {
    plot(assort_a, main = "Plot Network Property Assortativity", xlab = "Graph Index")
    abline(h=assort_truth, col = "red") 
    readline(prompt = "Pause. Press <Enter> to continue...")
  }
  
  print(length(which(assort_a >=  assort_truth))/length(assort_a))
  print(mean(assort_a))
  
  
  if (view_graphs == TRUE) {
    plot(density(assort_a), main = "Density plot Network Property Assortativity", xlab = "Assortativity Density")
    abline(v=assort_truth, col = "red") 
    readline(prompt = "Pause. Press <Enter> to continue...")
  }
  
  ####
  #x = sample(assort_a,10000, replace=TRUE, prob = assort_a^10)
  #mean(x)
  ####
  
  #	par(mfrow = c(ceiling(sqrt(max_deg)),ceiling(sqrt(max_deg))))
  #	for (i in c(1:(max_deg+1))) {
  #		plot(c(1:length(degree_dist_mat[,i])),cumsum(degree_dist_mat[,i]) / seq_along(degree_dist_mat[,i]), main = "Convergence of Degree Distribution", ylab = "Cum. Mean", xlab = "Iterations")
  #		abline(h=degree_dist_truth[i], col = "red") 
  #	}
  #	readline(prompt = "Pause. Press <Enter> to continue...")
  
  #	par(mfrow = c(ceiling(sqrt(max_deg)),ceiling(sqrt(max_deg))))
  #	for (i in c(1:(max_deg+1))) {
  #		plot(density(degree_dist_mat[,i]), main = "Density plot Network Property Degree Distribution", xlab = "Degree Density")
  #		abline(v=degree_dist_truth[i], col = "red") 
  #	}
  #	readline(prompt = "Pause. Press <Enter> to continue...")
  
  deg_dist_data = cbind(c(0:(max_deg)), degree_dist_truth, apply(degree_dist_mat,2,mean))
  colnames(deg_dist_data) = c("Degree", "TRUTH", "Estimated (Mean)")
  #	print(deg_dist_data)
  #	par(mfrow = c(1,1))
  boxplot(degree_dist_mat, main = "Degree Distribution", names = deg_dist_data[,1], xlab = "Node Degree")
  means= apply(degree_dist_mat,2,mean)
  means_true = degree_dist_truth
  points(1:length(means), means, pch = 23, cex = 0.75, bg = "blue")
  points(1:length(means_true), means_true, pch = 23, cex = 0.75, bg = "red")
  legend("topleft", legend=c("Estimated Degree Distribution (Mean)", "Truth"), pch=c(23,23), col = c("blue", "red")) 
  
  return(list(length(which(assort_a >=  assort_truth))/length(assort_a), mean(assort_a), assort_truth, mean(ecount_a), ecount_truth))
}