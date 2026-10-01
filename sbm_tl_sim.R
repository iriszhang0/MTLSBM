rm(list = ls())
source("/scratch/jz6169/MTLSBM/TransferLearning/sbm_tl_func.R")

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(task_id)

N <- 1000
###Vary target sparsity 
K <- 8; L <- 200; Kl_target <- 4; Kl_source <- sample(2:Kl_target, L, replace = T); rho_l_source <- rep(0.01,L)
rho_l_target_arr <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1)
results_sparsity <- vector("list", length(rho_l_target_arr))
for (i in seq_along(rho_l_target_arr)) {
  rho_l_target <- rho_l_target_arr[i]
  results_sparsity[[i]] <- run_simulation(N, K, L, Kl_target, Kl_source, rho_l_target, rho_l_source)
  results_sparsity[[i]]$rep <- task_id
}
results_sparsity <- data.table::rbindlist(results_sparsity)

###Vary number of layers
K <- 8; Kl_target <- 4; rho_l_target <- 0.05; L_arr <- c(25, 50, 75, 100, 125, 150, 175, 200)
results_layers <- vector("list", length(L_arr))
for (i in seq_along(L_arr)) {
  L <- L_arr[i]
  rho_l_source <- rep(0.01,L)
  # Kl_source <- rep(2,L); rho_l_source <- rep(0.01,L)
  Kl_source <- sample(2:Kl_target, L, replace = T) 
  results_layers[[i]] <- run_simulation(N, K, L, Kl_target, Kl_source, rho_l_target, rho_l_source)
  results_layers[[i]]$rep <- task_id
}
results_layers <- data.table::rbindlist(results_layers)


###Vary target community number
K <- 8; L <- 200; rho_l_source <- rep(0.01,L); rho_l_target <- 0.05; Kl_target_arr <- 2:K
results_communities <- vector("list", length(Kl_target_arr))
for (i in seq_along(Kl_target_arr)) {
  Kl_target <- Kl_target_arr[i]
  if(Kl_target == 2){
    Kl_source <- rep(2,L)
  }else{
    Kl_source <- sample(2:Kl_target, L, replace = T) 
  }
  results_communities[[i]] <- run_simulation(N, K, L, Kl_target, Kl_source, rho_l_target, rho_l_source)
  results_communities[[i]]$rep <- task_id
}
results_communities <- data.table::rbindlist(results_communities)

output_dir <- "/scratch/jz6169/MTLSBM/TransferLearning/Simulation/"

saveRDS(list(sparsity = results_sparsity, layers = results_layers, communities = results_communities),
        file = file.path(output_dir, paste0("result3_rep_", task_id, ".rds")))



