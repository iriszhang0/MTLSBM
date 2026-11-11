############################################################
# OpenFlights application (Section 5)
#
# Reproduces:
#   Figure 2  shared airport clusters (K = 4)
#   Figure 3  airline-specific airport clusters (K_l = 2)
#   Table 2   NLL of each airline's network under the 7 candidate
#             two-group splits of the shared clusters
#   Figure 4  route maps; routes are grouped by each airline's
#             lowest-NLL split from Table 2 (interpretation only)
#
# Nodes are the 100 busiest airports; each airline is one layer.
# Run from the folder containing routes.csv, airlines.csv and
# airports-extended.csv:   Rscript openflight_pipeline_clean.R
# Outputs are written to ./output/
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(mclust)          # adjustedRandIndex
  library(countrycode)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(patchwork)
  library(knitr)
})

set.seed(2026)
TOP_N <- 100   # number of airports
K     <- 4     # shared communities
KL    <- 2     # specific communities per airline

# Airlines reported in Table 2 and drawn in Figure 4
SELECTED_AIRLINES <- c(
  NH = "All Nippon Airways",
  AA = "American Airlines",
  DL = "Delta Air Lines",
  UA = "United Airlines",
  CA = "Air China",
  AF = "Air France"
)

OUT <- "output"
dir.create(OUT, showWarnings = FALSE)


############################################################
# 1. Data
############################################################

route_dat   <- fread("routes.csv", na.strings = "\\N")
airline_dat <- fread("airlines.csv", encoding = "UTF-8")
airport_dat <- fread("airports-extended.csv", encoding = "UTF-8")
colnames(airport_dat) <- c("id", "airport.name", "city.name", "country.name",
                           "IATA", "ICAO", "lat", "long", "altitude",
                           "tz.offset", "DST", "tz.name", "airport.type",
                           "source.data")

route_dat <- na.omit(route_dat)                      # drop incomplete records
route_dat[, `airline ID` := as.integer(`airline ID`)]
routes <- merge(route_dat, airline_dat[, .(`Airline ID`, Name)],
                by.x = "airline ID", by.y = "Airline ID", all.x = TRUE)
setnames(routes, c("source airport", "destination apirport"), c("src", "dst"))

cat(sprintf("Valid routes: %d | airlines: %d | airports: %d\n",
            nrow(route_dat), uniqueN(route_dat$`airline ID`),
            uniqueN(c(route_dat$`source airport`, route_dat$`destination apirport`))))

# Top-N airports by number of route endpoints (departures + arrivals)
top_airports <- routes[, .(airport = c(src, dst))][
  !is.na(airport), .N, by = airport][order(-N, airport)][seq_len(TOP_N), airport]
routes_top <- routes[src %in% top_airports & dst %in% top_airports &
                       src != dst & !is.na(Name)]

airport_nodes <- sort(unique(c(routes_top$src, routes_top$dst)))
airlines      <- sort(unique(routes_top$Name))
N <- length(airport_nodes)

# One symmetric adjacency matrix per airline
A_list <- lapply(setNames(airlines, airlines), function(a) {
  r <- routes_top[Name == a]
  A <- matrix(0, N, N, dimnames = list(airport_nodes, airport_nodes))
  idx <- cbind(match(r$src, airport_nodes), match(r$dst, airport_nodes))
  A[idx] <- 1
  A[idx[, 2:1, drop = FALSE]] <- 1
  A
})
cat(sprintf("Top-%d airports, %d airline layers\n", N, length(A_list)))

# Airport coordinates and geographic region (used to name clusters)
airport_info <- airport_dat[match(airport_nodes, IATA)]
stopifnot(!anyNA(airport_info$IATA))
airports <- data.table(Airport = airport_nodes,
                       airport.name = airport_info$airport.name,
                       country.name = airport_info$country.name,
                       lat = airport_info$lat, long = airport_info$long)
continent <- countrycode(airports$country.name, "country.name", "continent", warn = FALSE)
airports[, geo_region := fcase(country.name == "China", "China",
                               continent == "Americas", "Americas",
                               continent == "Europe", "Europe",
                               default = "Others")]


############################################################
# 2. Step 1: shared clusters (Figure 2)
############################################################

S0 <- matrix(0, N, N)
for (A in A_list) S0 <- S0 + A %*% A - diag(rowSums(A))

eg <- eigen(S0, symmetric = TRUE)
top <- order(eg$values, decreasing = TRUE)[1:K]
global_labels <- kmeans(eg$vectors[, top], centers = K, nstart = 50)$cluster

# Name each shared cluster by the region of most of its airports
cluster_name <- tapply(airports$geo_region, global_labels,
                       function(r) names(which.max(table(r))))
stopifnot(!anyDuplicated(cluster_name))
airports[, cluster := cluster_name[global_labels]]

cat("\nShared clusters vs geography (ARI =",
    round(adjustedRandIndex(airports$cluster, airports$geo_region), 3), ")\n")
print(airports[, table(cluster, geo_region)])

world <- ne_countries(scale = "medium", returnclass = "sf")
map_theme <- theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 16, hjust = 0.5),
        strip.text = element_text(size = 13))

fig2 <- ggplot(world) +
  geom_sf(fill = "gray95", color = "gray70") +
  geom_point(data = airports,
             aes(long, lat, color = factor(cluster,
                                           levels = c("Americas", "Europe", "China", "Others"),
                                           labels = c("Americas", "Europe", "Mainland China", "Others"))),
             size = 1.5) +
  labs(title = "Shared Airport Clusters", x = "Longitude", y = "Latitude", color = "Cluster") +
  map_theme
ggsave(file.path(OUT, "fig2_shared_clusters.pdf"), fig2, width = 10, height = 6)


############################################################
# 3. Step 2: airline-specific clusters (Figure 3)
############################################################

Theta <- matrix(0, N, K)
Theta[cbind(seq_len(N), global_labels)] <- 1
D <- diag(1 / colSums(Theta))

# Merge the K shared clusters into KL groups for one airline layer
specific_partition <- function(A) {
  B_hat <- D %*% t(Theta) %*% A %*% Theta %*% D
  eb <- eigen(B_hat)
  v <- eb$vectors[, order(abs(eb$values), decreasing = TRUE)[1:KL]]
  merged <- kmeans(v, centers = KL, nstart = 50)$cluster     # one label per shared cluster
  split(unname(cluster_name), merged)                        # regions in each group
}

# Readable name of a two-group split, e.g. "Europe vs Non-Europe"
split_name <- function(groups) {
  groups <- lapply(groups, function(g) sort(factor(g, levels = c("Americas", "Europe", "China", "Others"))))
  groups <- groups[order(lengths(groups), sapply(groups, function(g) as.integer(g[1])))]
  g1 <- paste(groups[[1]], collapse = "&"); g2 <- paste(groups[[2]], collapse = "&")
  if (length(groups[[1]]) == 1) paste0(g1, " vs Non-", g1) else paste(g1, "vs", g2)
}

specific <- rbindlist(lapply(names(A_list), function(a) {
  g <- specific_partition(A_list[[a]])
  g <- g[order(lengths(g))]
  data.table(airline = a, partition = split_name(g),
             group1 = paste(sort(g[[1]]), collapse = "&"))
}))
cat("\nSpecific partitions selected by Algorithm 1 (number of airlines):\n")
print(specific[, .N, by = partition][order(-N)])
fwrite(specific, file.path(OUT, "specific_partitions.csv"))

# One panel per distinct partition; Cluster 1 = the smaller group
panel_order <- c("Europe vs Non-Europe", "Americas vs Non-Americas",
                 "Others vs Non-Others", "China vs Non-China")
partitions <- unique(specific[, .(partition, group1)])
partitions <- partitions[order(match(partition, panel_order), partition)]
fig3_data <- rbindlist(lapply(seq_len(nrow(partitions)), function(i) {
  g1 <- strsplit(partitions$group1[i], "&")[[1]]
  airports[, .(long, lat, partition = partitions$partition[i],
               Cluster = factor(ifelse(cluster %in% g1, 1, 2)))]
}))
fig3_data[, partition := factor(partition, levels = partitions$partition)]

fig3 <- ggplot(world) +
  geom_sf(fill = "gray95", color = "gray70") +
  geom_point(data = fig3_data, aes(long, lat, color = Cluster), size = 0.6) +
  scale_color_manual(values = c("1" = "#F8766D", "2" = "#00BFC4")) +
  facet_wrap(~ partition, ncol = 2) +
  labs(title = "Airport Specific Clusters", x = "Longitude", y = "Latitude") +
  map_theme
ggsave(file.path(OUT, "fig3_specific_clusters.pdf"), fig3,
       width = 12, height = 3.2 * ceiling(nrow(partitions) / 2) + 1)


############################################################
# 4. Table 2: NLL of the 7 candidate splits
############################################################

# Group 1 = these clusters; group 2 = the remaining clusters
candidates <- list(
  "Americas"        = "Americas",
  "Europe"          = "Europe",
  "China"           = "China",
  "Others"          = "Others",
  "Americas&Europe" = c("Americas", "Europe"),
  "Americas&China"  = c("Americas", "China"),
  "Americas&Others" = c("Americas", "Others")
)

# Bernoulli NLL of an observed undirected network under a two-group SBM
# whose block probabilities are the observed densities (each pair once)
calculate_nll <- function(A, z, eps = 1e-15) {
  n1 <- sum(z == 1); n2 <- sum(z == 2)
  stopifnot(n1 >= 2, n2 >= 2)
  M <- matrix(c(choose(n1, 2), n1 * n2, n1 * n2, choose(n2, 2)), 2)
  C <- matrix(c(sum(A[z == 1, z == 1]) / 2, sum(A[z == 1, z == 2]),
                sum(A[z == 1, z == 2]),     sum(A[z == 2, z == 2]) / 2), 2)
  P <- pmin(pmax((C / M)[z, z], eps), 1 - eps)
  up <- upper.tri(A)
  -sum(A[up] * log(P[up]) + (1 - A[up]) * log1p(-P[up]))
}

stopifnot(all(SELECTED_AIRLINES %in% names(A_list)))
nll_table <- t(sapply(SELECTED_AIRLINES, function(a) sapply(candidates, function(cand)
  calculate_nll(A_list[[a]], ifelse(airports$cluster %in% cand, 1L, 2L)))))
rownames(nll_table) <- names(SELECTED_AIRLINES)

best <- data.table(
  code = rownames(nll_table),
  airline = unname(SELECTED_AIRLINES),
  best_candidate = colnames(nll_table)[apply(nll_table, 1, which.min)],
  best_nll = round(apply(nll_table, 1, min), 2))
best <- merge(best, specific[, .(airline, algorithm1_partition = partition)],
              by = "airline", sort = FALSE)
best[, best_partition := sapply(best_candidate, function(b)
  split_name(list(candidates[[b]], setdiff(cluster_name, candidates[[b]]))))]

cat("\nTable 2 (negative log-likelihood):\n")
print(round(nll_table, 2))
cat("\nLowest-NLL split vs Algorithm 1 partition:\n")
print(best[, .(code, best_partition, best_nll, algorithm1_partition)])

col_labels <- sapply(names(candidates), function(b)
  split_name(list(candidates[[b]], setdiff(cluster_name, candidates[[b]]))))
tab2 <- round(nll_table, 2); colnames(tab2) <- col_labels
fwrite(data.table(Airline = rownames(tab2), tab2), file.path(OUT, "table2_nll.csv"))
writeLines(kable(tab2, format = "latex", booktabs = TRUE, digits = 2,
                 caption = sprintf(paste(
                   "Negative log-likelihood of the observed network for seven candidate",
                   "specific community structures across %d selected airlines",
                   "(smaller values indicate better fit)."), nrow(tab2)),
                 label = "airline_loglikelihood"),
           file.path(OUT, "table2_nll.tex"))


############################################################
# 5. Figure 4: route maps, routes grouped by the lowest-NLL split
############################################################

region_fill <- c(Americas = "#F4ACB7", Europe = "#CAA8F5",
                 China = "#D4E09B", Others = "#A4B8C4")
route_cols <- c("#FE7F2D", "#00A6ED", "#FCCA46")   # within g1, between, within g2

route_map <- function(this_code) {
  airline <- SELECTED_AIRLINES[[this_code]]
  b <- best[code == this_code]
  g1 <- candidates[[b$best_candidate]]
  n1 <- paste(g1, collapse = "&")
  n2 <- if (length(g1) == 1) paste0("Non-", n1) else
    paste(intersect(c("Americas", "Europe", "China", "Others"),
                    setdiff(cluster_name, g1)), collapse = "&")

  A <- A_list[[airline]]
  e <- which(A == 1 & upper.tri(A), arr.ind = TRUE)
  in1 <- airports$cluster %in% g1
  levels_rt <- c(paste(n1, "-", n1), paste(n1, "-", n2), paste(n2, "-", n2))
  edges <- data.table(
    x = airports$long[e[, 1]], y = airports$lat[e[, 1]],
    xend = airports$long[e[, 2]], yend = airports$lat[e[, 2]],
    route = factor(fifelse(in1[e[, 1]] & in1[e[, 2]], levels_rt[1],
                   fifelse(!in1[e[, 1]] & !in1[e[, 2]], levels_rt[3], levels_rt[2])),
                   levels = levels_rt))
  nodes <- airports[sort(unique(c(e)))]

  ggplot(world) +
    geom_sf(fill = "gray80", color = "white") +
    geom_curve(data = edges, aes(x, y, xend = xend, yend = yend, color = route),
               curvature = 0.2, linewidth = 0.5) +
    geom_point(data = nodes, aes(long, lat, fill = cluster),
               shape = 21, size = 2, stroke = 0.2, color = "black") +
    scale_color_manual(values = setNames(route_cols, levels_rt), drop = FALSE,
                       name = "Route Type") +
    scale_fill_manual(values = region_fill, name = "Airport Cluster") +
    labs(title = sprintf("%s (%s): %s", airline, this_code, b$best_partition),
         x = "Longitude", y = "Latitude") +
    map_theme +
    theme(legend.box = "vertical", legend.title = element_text(size = 11),
          legend.text = element_text(size = 10)) +
    guides(fill = guide_legend(order = 1,
                               override.aes = list(shape = 21, size = 3, color = "black")),
           color = guide_legend(order = 2, override.aes = list(linewidth = 1.2)))
}

fig4 <- wrap_plots(lapply(names(SELECTED_AIRLINES), route_map), ncol = 2)
ggsave(file.path(OUT, "fig4_route_maps.pdf"), fig4,
       width = 18, height = 6 * ceiling(length(SELECTED_AIRLINES) / 2))

cat("\nDone. Outputs in", normalizePath(OUT), "\n")
