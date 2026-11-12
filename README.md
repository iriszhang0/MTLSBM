# MTLSBM: Community detection in sparse networks via transfer learning


## Repository contents

| File | Purpose |
|---|---|
| `sbm_tl_func.R` | Functions for Algorithm 1 and the simulation design: label generation, connection-probability matrices, the perturbation model, adjacency-matrix generation, spectral clustering, the target-layer merging step (`local_cluster`) and one simulation replicate (`run_simulation`). |
| `sbm_tl_sim.R` | Simulation study (Section 4): varies the target sparsity ρ_t, the number of layers L, the number of target communities K_t, and the perturbation proportion τ_ℓ. Runs one replicate per call; the replicate index is the random seed. |
| `openflight_pipeline_clean.R` | OpenFlights application (Section 5): preprocessing, shared airport clusters (Figure 2), airline-specific clusters (Figure 3), candidate-partition negative log-likelihoods (Table 2) and route maps (Figure 4). |
| `OpenFlightData/` | OpenFlights data: `routes.csv`, `airlines.csv`, `airports-extended.csv`, and `A_list.rds` (the 188 airline adjacency matrices on the top-100 airports). |

## Requirements

R (tested with 4.5.1) with the packages

```
Matrix, gtools, ggplot2, reshape2, mclust, aricode, dplyr, tidyr,
data.table, countrycode, sf, rnaturalearth, rnaturalearthdata,
patchwork, knitr
```

```r
install.packages(c("Matrix", "gtools", "ggplot2", "reshape2", "mclust",
                   "aricode", "dplyr", "tidyr", "data.table", "countrycode",
                   "sf", "rnaturalearth", "rnaturalearthdata", "patchwork",
                   "knitr"))
```

## OpenFlights application (Figures 2–4, Table 2)

`openflight_pipeline_clean.R` reads the three CSV files from the working
directory, so run it from `OpenFlightData/`:

```sh
cd OpenFlightData
Rscript ../openflight_pipeline_clean.R
```

Outputs are written to `OpenFlightData/output/`:

| Result | Output file |
|---|---|
| Figure 2: shared airport clusters (K = 4) | `fig2_shared_clusters.pdf` |
| Figure 3: airline-specific airport clusters (K_ℓ = 2) | `fig3_specific_clusters.pdf` |
| Table 2: negative log-likelihood of each candidate partition | `table2_nll.csv`, `table2_nll.tex` |
| Figure 4: route maps of the selected airlines | `fig4_route_maps.pdf` |
| Algorithm 1 partition for every airline | `specific_partitions.csv` |

What the script does:

1. **Preprocessing.** Drops incomplete route records (66,765 routes, 547
   airlines and 3,246 airports remain), keeps routes between the 100
   airports with the most route endpoints, and builds one symmetric
   adjacency matrix per airline (188 layers). It prints these counts.
2. **Shared clusters (Step 1 of Algorithm 1).** Spectral clustering of the
   bias-adjusted matrix S = Σ_ℓ (A_ℓ² − D_ℓ) with K = 4. Each cluster is
   named after the region most of its airports belong to (Americas,
   Europe, mainland China, Others); the script prints the agreement with
   geography (ARI).
3. **Specific clusters (Step 2).** Each airline in turn is the target
   layer; its 4 shared clusters are merged into K_ℓ = 2 groups using the
   leading eigenvectors of B̂_ℓ.
4. **Table 2.** For each selected airline and each of the 7 ways to split
   the 4 shared clusters into two groups, the block probabilities are the
   observed route densities, and the Bernoulli negative log-likelihood of
   the airline's observed network is computed over all 4,950 airport pairs
   (each pair counted once).
5. **Figure 4.** Each selected airline's routes are grouped by its
   lowest negative log-likelihood split from Table 2. The figure is meant
   for interpretation.

The airlines in Table 2 and Figure 4 are set by `SELECTED_AIRLINES` at the
top of the script (NH, AA, DL, UA, CA, AF). The random seed is fixed with
`set.seed(2026)`.

## Simulation study (Section 4)

`sbm_tl_sim.R` runs one replicate (n = 1,000 nodes) of four experiments:

| Experiment | Settings | Output element |
|---|---|---|
| Target sparsity ρ_t (Figure 1, upper left) | ρ_t = 0.01, …, 0.10; K = 8, L = 200, K_t = 4, ρ_ℓ = 0.01 | `sparsity` |
| Number of layers L (Figure 1, upper right) | L = 25, …, 200; K = 8, K_t = 4, ρ_t = 0.05 | `layers` |
| Target communities K_t (Figure 1, lower left) | K_t = 2, …, 8; K = 8, L = 200, ρ_t = 0.05 | `communities` |
| Perturbation τ_ℓ (Figure 1, lower right) | τ_ℓ = 0, 0.02, …, 0.20; K = 8, K_t = 4, L = 200, ρ_t = 0.05 | `permutation` |

It was run as a SLURM array of 100 tasks (`SLURM_ARRAY_TASK_ID` = 1, …, 100
is the replicate index and the seed), with R 4.5.1 and one CPU per task.
Each `run_simulation()` call takes on the order of 10–20 seconds on one
core.

## K and K_t estimation 

Estimation of K and K_t with the Python package `kneed` with curve="convex" and direction="decreasing" on the sorted absolute value of eigenvalues. 


## Data

The OpenFlights data are from https://www.kaggle.com/datasets/open-flights
(original source: https://openflights.org/data, Open Database License).
