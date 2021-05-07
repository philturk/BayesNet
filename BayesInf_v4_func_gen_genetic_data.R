########Phylogenetic Data Generate############33

Genetic_Seq_Data <- function(P, Ia, final_time = max(Ia[which(Ia < Inf)]), genetic_bits = 128) {
  
  infection_order = order(Ia)
  Genetic_Seq_mat = matrix(0, nrow = length(Ia), ncol = genetic_bits)
  update_times_a = array(0, dim = length(Ia))
  Genetic_Seq_mat = cbind(Genetic_Seq_mat, update_times_a)
  for(i in infection_order[c(1:length(which(Ia < Inf)))]) {
    parent_i = which(P[,i]==1)
    if (length(parent_i) == 0) { #root node
      Genetic_Seq_mat[i,] = c(genetic_seq_root_node(genetic_bits), Ia[i])
    } else { #has a parent
      Genetic_Seq_mat[parent_i,] = c(genetic_seq_parent_node(Genetic_Seq_mat[parent_i,-(genetic_bits+1)], Genetic_Seq_mat[parent_i,(genetic_bits+1)], Ia[i], genetic_bits), Ia[i])
      Genetic_Seq_mat[i,] = c(genetic_seq_inf_node(Genetic_Seq_mat[parent_i,-(genetic_bits+1)], Genetic_Seq_mat[parent_i,(genetic_bits+1)], Ia[i], genetic_bits), Ia[i])
    }
  }
  
  ###Update everyone to be sampled at end of trial
  
  for(i in infection_order[c(1:length(which(Ia < Inf)))]) {
    Genetic_Seq_mat[i,] = c(genetic_seq_parent_node(Genetic_Seq_mat[i,-(genetic_bits+1)], Genetic_Seq_mat[i,(genetic_bits+1)], final_time, genetic_bits), final_time)
  }
  return(Genetic_Seq_mat)
  
  
}


genetic_seq_root_node <- function(genetic_bits) {
  return(sample(c(0,1), genetic_bits, replace = TRUE))
}

genetic_seq_parent_node <- function(Genetic_Seq_a, last_updatetime, ctime, genetic_bits) {
  nbit = round(ctime - last_updatetime)
  change_bits = array(0, dim = genetic_bits)
  change_bits[sample(c(1:genetic_bits), nbit, replace = FALSE)] = 1
  Genetic_Seq_a = (Genetic_Seq_a + change_bits) %% 2
  return(Genetic_Seq_a)
}

genetic_seq_inf_node <- function(Genetic_Seq_a, last_updatetime, ctime, genetic_bits) {
  nbit = 0
  change_bits = array(0, dim = genetic_bits)
  change_bits[sample(c(1:genetic_bits), nbit, replace = FALSE)] = 1
  Genetic_Seq_a = (Genetic_Seq_a + change_bits) %% 2
  return(Genetic_Seq_a)
}