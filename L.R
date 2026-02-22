rm(list = ls())
library(Matrix)
library(gtools) 
library(ggplot2)
library(reshape2)
library(mclust)
library(aricode)

############Functions###############
## Generate True Global Membership 
gen_true_labels <- function(K,N){
  groups <- rep(1:K, each = ceiling(N / K), length.out = N)
  return(sort(groups))
}

## Generate True Local Membership 
gen_true_label_list <- function(K,N,Kl,L,true_labels){
  base_labels <- true_labels
  label_list <- vector("list", L)
  base_size <- K %/% Kl
  remainder <- K %% Kl
  group_sizes_template <- rep(base_size, Kl)
  if (remainder > 0) {
    group_sizes_template[1:remainder] <- 
      group_sizes_template[1:remainder] + 1
  }
  for (i in seq_len(L)) {
    permuted_classes <- sample(1:K)
    group_sizes <- sample(group_sizes_template)
    class_assignment <- numeric(K)
    start_idx <- 1
    for (g in 1:Kl) {
      end_idx <- start_idx + group_sizes[g] - 1
      class_assignment[permuted_classes[start_idx:end_idx]] <- g
      start_idx <- end_idx + 1
    }
    new_labels <- class_assignment[base_labels]
    
    label_list[[i]] <- new_labels
  }
  return(label_list)
}

## Generate Probability Matrix For each Layer 
B_matrix_func <- function(Kl, rho, true_label_list, n_rows = N){
  B_mat <- list()
  for (t in 1:length(true_label_list)){
    true_label <- true_label_list[[t]]
    if (t%%2==0){
      block_entry <- c(0.8, 0.8, 0.3)
    }else{
      block_entry <- c(0.3, 0.3, 0.8)
    }
    B <- matrix(NA,nrow=n_rows,ncol=n_rows)
    for (l in 1:Kl) {
      ind <- which(true_label == l) 
      B[ind, ind] <- block_entry[l] 
    }
    combs <- combn(Kl, 2) 
    counter <- Kl + 1  
    for (k in 1:ncol(combs)) {
      i <- combs[1,k]
      j <- combs[2,k]
      value <- block_entry[counter]
      rows_i <- which(true_label == i)
      rows_j <- which(true_label == j)
      B[rows_i, rows_j] <- value
      B[rows_j, rows_i] <- value
      counter <- counter + 1
    }
    B_mat[[t]] <- B*rho
  }
  return (B_mat)
}

## Shared Clustering (SoS & SoSDebiased)
global_cluster_func <- function(matrix_list, true_labels, K) {
  clustering_results <- list()
  adj_rand_ind <- numeric(length(matrix_list))
  nor_mut_info <- numeric(length(matrix_list))
  for (matrix_idx in seq_along(matrix_list)) {
    rho_key <- names(matrix_list)[matrix_idx]
    mat <- matrix_list[[rho_key]]
    eigen_decomp <- eigen(mat)
    eigen_values <- eigen_decomp$values
    eigen_top_loc <- order(eigen_values, decreasing = TRUE)[1:K]
    top_eigenvectors <- eigen_decomp$vectors[, eigen_top_loc]
    clustering <- kmeans(top_eigenvectors, centers = K, nstart = 10)
    predicted_labels <- clustering$cluster
    clustering_results[[rho_key]] <- predicted_labels
    ari <- adjustedRandIndex(true_labels, predicted_labels)
    nmi <- NMI(true_labels, predicted_labels)
    adj_rand_ind[matrix_idx] <- ari
    nor_mut_info[matrix_idx] <- nmi
  }
  list(clustering_results = clustering_results, 
       adj_rand_ind = adj_rand_ind,
       nor_mut_info = nor_mut_info) 
}


## Specific Clustering (SoS & SoSDebiased)
local_cluster_func <- function(A_mat, pred_labels, true_label_list, method, K) {
  adj_rand_ind_layers <- list()
  nor_mut_info_layers <- list()
  avg_adj_rand_ind_local <- numeric(length(A_mat))
  avg_nor_mut_info_local <- numeric(length(A_mat))
  for (prop_idx in seq_along(A_mat)) {
    prop_key <- names(A_mat)[prop_idx]
    adj_rand_ind_layers[[prop_key]] <- list()
    nor_mut_info_layers[[prop_key]] <- list()
    A_l_list <- A_mat[[prop_idx]]  
    theta_hat <- pred_labels[[prop_idx]]
    one_hot <- matrix(0, nrow = length(theta_hat), ncol = K)
    one_hot[cbind(seq_along(theta_hat), theta_hat)] <- 1
    theta_hat <- one_hot
    D <- diag(rep(0,K))
    for (i in 1:K){
      D[i,i] <- 1/sum(pred_labels[[prop_idx]]==i)
    }
    for (layer_idx in seq_along(A_mat[[prop_idx]])){
      A_l <- A_l_list[[layer_idx]]
      true_labels <- true_label_list[[layer_idx]]
      layer_dat <- D %*% t(theta_hat) %*% A_l %*% theta_hat %*% D
      cluster_num <- length(unique(true_labels))
      if(cluster_num == nrow(layer_dat)){
        new_pred_labels <- pred_labels[[prop_idx]] 
      }else{
        if(method == "kmeans"){
          clustering <- kmeans(layer_dat, centers = cluster_num, nstart = 2)
        }
        if(method == "spectral"){
          eigen_decomp <- eigen(layer_dat)
          eigen_values <- eigen_decomp$values
          eigen_top_loc <- order(abs(eigen_values), decreasing = TRUE)[1:cluster_num]
          top_eigenvectors <- eigen_decomp$vectors[, eigen_top_loc]
          clustering <- kmeans(top_eigenvectors, centers = cluster_num, nstart = 2)
        }
        predicted_labels <- clustering$cluster
        pred_labels_og <- pred_labels[[prop_idx]]
        unique_pred_labels <- 1:K 
        unique_predicted_labels <- predicted_labels
        transformation_rules <- setNames(unique_predicted_labels, as.character(unique_pred_labels))
        new_pred_labels <- as.integer(sapply(pred_labels_og, function(x) transformation_rules[as.character(x)]))
      }
      ari <- adjustedRandIndex(true_labels, new_pred_labels)
      nmi <- NMI(true_labels, new_pred_labels)
      adj_rand_ind_layers[[prop_key]][[layer_idx]] <- ari
      nor_mut_info_layers[[prop_key]][[layer_idx]] <- nmi
    }
    avg_adj_rand_ind_local[prop_idx] <- mean(unlist(adj_rand_ind_layers[[prop_key]]))
    avg_nor_mut_info_local[prop_idx] <- mean(unlist(nor_mut_info_layers[[prop_key]]))
  }
  list(
    avg_adj_rand_ind_local = avg_adj_rand_ind_local, 
    avg_nor_mut_info_local = avg_nor_mut_info_local) 
}

## Specific Clustering (Single Layer)
single_layer_func <- function(A_list, Kl, true_label_list){
  L = length(A_list)
  adj_rand_ind_local <- numeric(L)
  for (l in 1:L){
    A = A_list[[l]]
    true_labels <- true_label_list[[l]]
    eigen_decomp = eigen(A)
    eigen_values <- eigen_decomp$values
    eigen_top_loc <- order(abs(eigen_values), decreasing = TRUE)[1:Kl]
    top_eigenvectors <- eigen_decomp$vectors[, eigen_top_loc]
    clustering <- kmeans(top_eigenvectors, centers = Kl, nstart = 50)
    predicted_labels <- clustering$cluster
    adj_rand_ind_local[l] <- adjustedRandIndex(true_labels, predicted_labels)
  }
  return(mean(adj_rand_ind_local))
}


run_simulation <- function(seed, K, Kl, rho, prop, L_values) {
  set.seed(seed)
  results <- data.frame(
    K = numeric(0),
    Kl = numeric(0), 
    rho = numeric(0), 
    prop = numeric(0), 
    L = numeric(0), 
    S0_ari_global = numeric(0),
    M_ari_global = numeric(0),
    S0_ari_local_spectral = numeric(0),
    M_ari_local_spectral = numeric(0),
    local_ari_single_layer = numeric(0)
  )
  for (L in L_values) {
    true_labels <- gen_true_labels(K, N) 
    true_label_list <- gen_true_label_list(K, N, Kl, L, true_labels) 
    B_mat <- B_matrix_func(Kl, rho, true_label_list, N)
    A_list <- lapply(B_mat, function(B) {
      A <- matrix(0, nrow = nrow(B), ncol = ncol(B))
      upper_tri_indices <- which(upper.tri(A))
      A[upper_tri_indices] <- rbinom(length(upper_tri_indices), size = 1, prob = B[upper_tri_indices])
      diag(A) <- 0
      A_symmetric <- A + t(A)
      return(A_symmetric)
    })
    local_ari_single_layer = single_layer_func(A_list, Kl, true_label_list)
    S0_mat <- list()
    M_mat <- list()
    n <- nrow(A_list[[1]])
    S0 <- matrix(0, nrow = n, ncol = n)
    M <- matrix(0, nrow = n, ncol = n)
    for (A in A_list) {
      A_squared <- A %*% A
      degrees <- rowSums(A)
      D <- diag(degrees) 
      S0 <- S0 + (A_squared - D)
      M <- M + A_squared 
    }
    S0_mat[["Prop_0"]] <- S0
    M_mat[["Prop_0"]] <- M
    S0_results <- global_cluster_func(S0_mat, true_labels, K)
    M_results <- global_cluster_func(M_mat, true_labels, K)
    S0_labels <- S0_results$clustering_results
    M_labels <- M_results$clustering_results
    A_mat <- list()
    A_mat[["Prop_0"]]  <- A_list
    S0_layers_spectral <- local_cluster_func(A_mat, S0_labels, true_label_list, method = "spectral", K)
    M_layers_spectral <- local_cluster_func(A_mat, M_labels, true_label_list, method = "spectral", K)
    
    results <- rbind(results, data.frame(
      K = K,
      Kl = Kl, 
      rho = rho, 
      prop = prop, 
      L = L, 
      S0_ari_global = S0_results$adj_rand_ind,
      M_ari_global = M_results$adj_rand_ind,
      S0_ari_local_spectral = S0_layers_spectral$avg_adj_rand_ind_local,
      M_ari_local_spectral = M_layers_spectral$avg_adj_rand_ind_local,
      local_ari_single_layer = local_ari_single_layer
    ))
  }
  return(results)
}




###########Parameters############
K <- 16
Kl <- 2
N <- 1000
rho <- 0.02
permu_prop <- 0
L_values <- c(50,100,150,200,250,300,350,400)

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
result <- run_simulation(seed = task_id, 
                         K = K, 
                         Kl = Kl, 
                         rho = rho, 
                         prop=permu_prop, 
                         L_values = L_values)

output_dir <- "/scratch/jz6169/MTLSBM/L/"

result_file <- paste0(output_dir, "result_df_", task_id, ".rds")
saveRDS(result, file = result_file)


