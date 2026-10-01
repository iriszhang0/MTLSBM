###########################################
source("/scratch/jz6169/MTLSBM/TransferLearning/sbm_tl_func.R")
library(Matrix)
library(ggplot2)
library(reshape2)
library(aricode) 
library(dplyr)
library(tidyr)
library(data.table)


compute_test_nll <- function(A_train, A_test, labels) {
  N <- nrow(A_train)
  K <- max(labels)
  P_hat <- matrix(0, K, K)
  for (a in 1:K) {
    for (b in a:K) {
      idx_a <- which(labels == a)
      idx_b <- which(labels == b)
      if (a == b) {
        # Within-community pairs
        pairs <- combn(idx_a, 2)
        edges <- sum(A_train[pairs[1, ], pairs[2, ]])
        n_pairs <- ncol(pairs)
      } else {
        # Between-community pairs
        edges <- sum(A_train[idx_a, idx_b])
        n_pairs <- length(idx_a) * length(idx_b)
      }
      # MLE of block probability
      P_hat[a, b] <- edges / n_pairs
      P_hat[b, a] <- P_hat[a, b]
    }
  }
  # Avoid log(0)
  eps <- 1e-10
  P_hat <- pmin(pmax(P_hat, eps), 1 - eps)
  # Compute test NLL using only i < j
  nll <- 0
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      p_ij <- P_hat[labels[i], labels[j]]
      y_ij <- A_test[i, j]
      nll <- nll - (y_ij * log(p_ij) + (1 - y_ij) * log(1 - p_ij))
    }
  }
  return(nll)
}


run_layer_sampling <- function(A_list, L_target, K_est = 4, Kl_est_value = 2, train_prop = 0.9) {
  L_all <- length(A_list)
  selected_layers <- sample(seq_len(L_all), L_target, replace = FALSE)
  A_selected <- A_list[selected_layers]
  L <- length(A_selected)
  N <- nrow(A_selected[[1]])
  Kl_est <- rep(Kl_est_value, L)
  A_list_train <- vector("list", L)
  A_list_test <- vector("list", L)
  idx <- which(upper.tri(matrix(0, N, N)), arr.ind = TRUE)
  for (l in seq_len(L)) {
    T <- matrix(0, N, N)
    sampled <- runif(nrow(idx)) < train_prop
    T[cbind(idx[, 1], idx[, 2])] <- sampled
    T[cbind(idx[, 2], idx[, 1])] <- sampled
    A_list_train[[l]] <- A_selected[[l]] * T
    A_list_test[[l]] <- A_selected[[l]] * (1 - T)
  }
  S0_train <- matrix(0, N, N)
  S0_test <- matrix(0, N, N)
  for (l in seq_len(L)) {
    A <- A_list_train[[l]]
    S0_train <- S0_train + A %*% A - diag(rowSums(A))
    A <- A_list_test[[l]]
    S0_test <- S0_test + A %*% A - diag(rowSums(A))
  }
  global_clustering_train <- spectral_clustering(S0_train, K_est, true_labels = NA)
  pred_global_labels_train <- global_clustering_train$predicted_labels
  pred_local_labels_train <- vector("list", L)
  target_only_labels_train <- vector("list", L)
  for (l in seq_len(L)) {
    current_list <- c(A_list_train[l:L], A_list_train[1:(l - 1)])
    pred_local_labels_train[[l]] <- local_cluster(current_list, pred_global_labels_train, 
                                                  true_target_labels = NA, K_est, Kl_est[l])$predicted_labels
    target_only_labels_train[[l]] <- spectral_clustering(current_list[[1]], Kl_est[l], true_labels = NA)$predicted_labels
  }
  test_nll_local <- numeric(L)
  test_nll_target_only <- numeric(L)
  for (l in seq_len(L)) {
    test_nll_local[l] <- compute_test_nll(A_train = A_list_train[[l]],
                                          A_test = A_list_test[[l]],
                                          labels = pred_local_labels_train[[l]])
    test_nll_target_only[l] <- compute_test_nll(A_train = A_list_train[[l]], 
                                                A_test = A_list_test[[l]], 
                                                labels = target_only_labels_train[[l]])
  }
  return(list(L = L_target, mean_nll_local = mean(test_nll_local), 
              mean_nll_target_only = mean(test_nll_target_only)))
}


task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(task_id)
S0 <- readRDS("/scratch/jz6169/MTLSBM/TransferLearning/OpenFlight/S0_mat.rds")
A_list <- readRDS("/scratch/jz6169/MTLSBM/TransferLearning/OpenFlight/A_list.rds")
L_values <- c(25,50,75,100,125,150,188)
results <- lapply(L_values, function(L_target) {
  run_layer_sampling(A_list = A_list, L_target = L_target)
  })


output_dir <- "/scratch/jz6169/MTLSBM/TransferLearning/OpenFlight/"

saveRDS(results, 
        file = file.path(output_dir, paste0("result_L_", task_id, ".rds")))




