rm(list = ls())
library(Matrix)
library(gtools) 
library(ggplot2)
library(reshape2)
library(mclust)
library(aricode) 
library(dplyr)
library(tidyr)
library(data.table)
# library(kneedle)

gen_true_labels <- function(K,N){
  groups <- rep(1:K, each = ceiling(N / K), length.out = N)
  return(sort(groups))
}


gen_true_label_list <- function(K,N,Kl_source,L,true_labels){
  base_labels <- true_labels
  label_list <- vector("list", L)
  for (i in seq_len(L)) {
    Kl <- Kl_source[i]
    base_size <- K %/% Kl
    remainder <- K %% Kl
    group_sizes_template <- rep(base_size, Kl)
    if (remainder > 0) {
      group_sizes_template[1:remainder] <- group_sizes_template[1:remainder] + 1
    }
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


B_matrix_func <- function(Kl_source, rho_source, true_local_label, n_rows = N) {
  B_mat <- list()
  for (t in seq_along(true_local_label)) {
    true_label <- true_local_label[[t]]
    Kl <- Kl_source[t]
    rho <- rho_source[t]
    if (t %% 2 == 1) {
      within_prob <- 0.8
      cross_prob  <- 0.3
    } else {
      within_prob <- 0.3
      cross_prob  <- 0.8
    }
    B_block <- matrix(cross_prob, nrow = Kl, ncol = Kl)
    diag(B_block) <- within_prob
    B <- B_block[true_label, true_label]
    B_mat[[t]] <- B * rho
  }
  return(B_mat)
}


B_matrix_permutation <- function(B_mat_source, true_source_labels, Kl_source, rho_l_source, permu_prop) {
  layers_list <- list()
  new_labels <- true_source_labels
  permuted_count <- integer(N)
  for (t in seq_along(B_mat_source)) {
    labels <- true_source_labels[[t]]
    N <- length(labels)
    Kl <- Kl_source[t]
    num_changes <- ceiling(permu_prop * N)
    nodes_to_change <- sample(seq_len(N), num_changes)
    permuted_count[nodes_to_change] <- permuted_count[nodes_to_change] + 1
    for (node in nodes_to_change) {
      old_group <- labels[node]
      new_group <- sample(setdiff(1:Kl, old_group), 1)
      new_labels[[t]][node] <- new_group
    }
  }
  new_layers <- B_matrix_func(Kl_source, rho_l_source, new_labels)
  layers_list <- new_layers
  max_permutation_prop <- max(permuted_count / L)
  list(B_mat_source = layers_list, kappa_max = max_permutation_prop)
}




A_mat_func <- function(B_mat) {
  lapply(B_mat, function(B) {
    A <- matrix(0, nrow = nrow(B), ncol = ncol(B))
    upper_tri_indices <- which(upper.tri(A))
    A[upper_tri_indices] <- rbinom(length(upper_tri_indices), size = 1, prob = B[upper_tri_indices])
    diag(A) <- 0
    A_symmetric <- A + t(A)
    return(A_symmetric)
  })
}



spectral_clustering <- function(adj_matrix, cluster_num, true_labels){
  eigen_decomp <- eigen(adj_matrix)
  eigen_values <- eigen_decomp$values
  eigen_top_loc <- order(abs(eigen_values), decreasing = TRUE)[1:cluster_num]
  top_eigenvectors <- eigen_decomp$vectors[, eigen_top_loc]
  clustering <- kmeans(top_eigenvectors, centers = cluster_num, nstart = 10)
  predicted_labels <- clustering$cluster
  if(length(true_labels) == 1 && is.na(true_labels)){
    ari = NA
  }else{
    ari <- adjustedRandIndex(true_labels, predicted_labels)
  }
  list(ari = ari, predicted_labels = predicted_labels)
}


local_cluster <- function(A_list, pred_global_labels, true_target_labels, K, Kl_target) {
  theta_hat <- pred_global_labels
  one_hot <- matrix(0, nrow = length(theta_hat), ncol = K)
  one_hot[cbind(seq_along(theta_hat), theta_hat)] <- 1
  theta_hat <- one_hot
  D <- diag(rep(0,K))
  for (i in 1:K){
    D[i,i] <- 1/sum(pred_global_labels==i)
  }
  A_l <- A_list[[1]]
  layer_dat <- D %*% t(theta_hat) %*% A_l %*% theta_hat %*% D
  cluster_num <- Kl_target
  if(cluster_num == nrow(layer_dat)){
    new_pred_labels <- pred_global_labels
  }else{
    eigen_decomp <- eigen(layer_dat)
    eigen_values <- eigen_decomp$values
    eigen_top_loc <- order(abs(eigen_values), decreasing = TRUE)[1:cluster_num]
    top_eigenvectors <- eigen_decomp$vectors[, eigen_top_loc]
    clustering <- kmeans(top_eigenvectors, centers = cluster_num, nstart = 2)
    predicted_labels <- clustering$cluster
    pred_labels_og <- pred_global_labels
    unique_pred_labels <- 1:K 
    unique_predicted_labels <- predicted_labels
    transformation_rules <- setNames(unique_predicted_labels, as.character(unique_pred_labels))
    new_pred_labels <- as.integer(sapply(pred_labels_og, function(x) transformation_rules[as.character(x)]))
  }
  if(length(true_target_labels) == 1 && is.na(true_target_labels)){
    ari <- NA
  }else{
    true_labels <- true_target_labels[[1]]
    ari <- adjustedRandIndex(true_labels, new_pred_labels)
  }
  list(ari = ari, predicted_labels = new_pred_labels)
}


shared_community_num_estimate <- function(A_list, S0, N, L, prop = 0){
  K_estimate <- NA
  sum_A_mat <- sum(sapply(A_list,sum)) 
  rho_hat <- sum_A_mat/((L)*N*N)
  if(prop == 0){
      t_n <- (L)*N*rho_hat^2+sqrt(L)*N*rho_hat*(log(L+N))^(0.5)
  }else{
      t_n <- L*N*rho_hat^2*sqrt(log(N))+sqrt(L)*N*rho_hat*log(L+N)+L*N^2*rho_hat^2*(prop)^(5/6)
  }
  eigen_decomp <- eigen(S0)
  eigen_values <- eigen_decomp$values
  eigen_values_top <- sort(eigen_values, decreasing = TRUE)
  K_estimate <- which.max(eigen_values_top>t_n & c(eigen_values_top[-1],-Inf)<t_n)
  return(K_estimate)
}


target_community_num_estimate <- function(A_list, predicted_labels, K_est, prop = 0){
  Kt_estimate <- NA
  A_t <- A_list[[1]]
  theta_hat <- predicted_labels
  one_hot <- matrix(0, nrow = length(theta_hat), ncol = K_est)
  one_hot[cbind(seq_along(theta_hat), theta_hat)] <- 1
  theta_hat <- one_hot
  D <- diag(rep(0,K_est))
  for (i in 1:K_est){
      D[i,i] <- 1/sum(predicted_labels==i)
  }
  layer_dat <- D %*% t(theta_hat) %*% A_t %*% theta_hat %*% D
  rho_hat_t <- sum(A_t) / (N*N)
  if(prop == 0){
    t_t_n <- K_est^(2.1)*sqrt(rho_hat_t*log(N))/N
  }else{
    t_t_n <- K_est^(3/2)*rho_hat_t*prop+K_est^(5/2)*sqrt(rho_hat_t*log(N))/N
  }
  eigen_decomp <- eigen(layer_dat)
  eigen_values <- eigen_decomp$values
  eigen_values_top <- sort(abs(eigen_values), decreasing = TRUE)
  Kt_estimate <- which.max(eigen_values_top>t_t_n & c(eigen_values_top[-1],-Inf)<t_t_n)
  return(Kt_estimate)
}



run_simulation <- function(N, K, L, Kl_target, Kl_source, rho_l_target, rho_l_source, permu_prop = 0){
  true_global_labels <- gen_true_labels(K, N) 
  true_target_labels <- gen_true_label_list(K, N, Kl_target, 1, true_global_labels)
  true_source_labels <- gen_true_label_list(K, N, Kl_source, L, true_global_labels)
  B_mat_target <- B_matrix_func(Kl_target, rho_l_target, true_target_labels, n_rows = N)
  B_mat_source <- B_matrix_func(Kl_source, rho_l_source, true_source_labels, n_rows = N)
  kappa_max <- 0 
  if(permu_prop != 0){
    B_mat_pertub_func <- B_matrix_permutation(B_mat_source, true_source_labels, Kl_source, rho_l_source, permu_prop)
    B_mat_source <- B_mat_pertub_func$B_mat_source
    kappa_max <- B_mat_pertub_func$kappa_max
  }
  A_target <- A_mat_func(B_mat_target)
  A_source <- A_mat_func(B_mat_source)
  #single layer 
  target_only_ari <- spectral_clustering(A_target[[1]], Kl_target, true_target_labels[[1]])$ari
  #proposed method 
  S0 <- matrix(0, nrow = N, ncol = N)
  A_list <- vector("list", L + 1)
  A_list[[1]] <- A_target[[1]]
  A_list[2:(L+1)] <- A_source 
  for (A in A_list) {
    A_squared <- A %*% A
    degrees <- rowSums(A)
    D <- diag(degrees) 
    S0 <- S0 + (A_squared - D)
  }
  K_estimated <- shared_community_num_estimate(A_list, S0, N, L, permu_prop)
  global_clustering <- spectral_clustering(S0, K, true_global_labels)
  pred_global_labels <- global_clustering$predicted_labels
  global_ari <- global_clustering$ari
  K_estimate <- max(1, min(K_estimate, K))
  Kl_target_estimated <- target_community_num_estimate(A_list, pred_global_labels, K_estimated, permu_prop)
  target_ari <- local_cluster(A_list, pred_global_labels, true_target_labels, K, Kl_target)$ari
  list(K = K, L = L, Kl_target = Kl_target, rho_l_target = rho_l_target, permu_prop = permu_prop, 
    target_only_ari = target_only_ari, 
    global_sbm_ari = global_ari,
    target_sbm_ari = target_ari, 
    K_estimated = K_estimated, 
    Kl_target_estimated = Kl_target_estimated, 
    kappa_max = kappa_max)
}








