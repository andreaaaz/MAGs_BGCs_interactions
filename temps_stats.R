################################################
### Network stats by temperature  ######################
## Andrea Zermeño Díaz #########################
# october-2026 ###############################

# compare the effect of the temperature on the networks

# libraries
library(igraph)
library(tidyverse)
library(ggplot2)

# data load 
setwd("interactions/2026-09-interactions08/")
temperature <- c("global", "high", "mid", "low")
types <- c("MAG-BGC", "MAG-MAG-rec", "MAG-MAG")

make_network <- function(temperature, type) {
  if (type == "MAG-BGC") {
    folder <- "mOTUs_Species_Cluster_gcc"
    file_type <- "mb"
  } else if (type == "MAG-MAG-rec") {
    folder <- "mOTUs_Species_Cluster_gcc"
    file_type <- "mm"
  } else if (type == "MAG-MAG") {
    folder <- "mOTUs_Species_Cluster"
    file_type <- "mm"
  } 

  edges <- read.csv(file.path(folder, temperature, paste0("edges_", file_type, ".csv")))
  nodes <- read.csv(file.path(folder, temperature, paste0("nodes_", file_type, ".csv")))
  graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
}

networks_mb <- list()
networks_mmr <- list()
networks_mm <- list()

for (temperature in temperature) {
  networks_mb[[temperature]] <- make_network(temperature, "MAG-BGC")
  networks_mmr[[temperature]] <- make_network(temperature, "MAG-MAG-rec")
  networks_mm[[temperature]] <- make_network(temperature, "MAG-MAG")
}

# NODE STATS


# calculate degree, betweeennes, closeness and eigenvector
get_node_stats <- function(g, network_name, network_type) {
  tibble(network = network_name,
         network_type = network_type,
         temperature = network_name,
         node = V(g)$name,
         degree = degree(g),
         norm_degree = degree(g, normalized = TRUE),
         betweenness = betweenness(g, weights = NA),
         norm_betweenness = betweenness(g, weights = NA, normalized = TRUE),
         closeness = closeness(g, weights = NA),
         norm_closeness = closeness(g, weights = NA, normalized = TRUE),
         eigenvector = eigen_centrality(g, weights = NA)$vector)
}

stats_mb <- purrr::imap_dfr(networks_mmr, ~ get_node_stats(.x, .y, "MAG-MAG-rec"))

densityplot <- function(node_stats, stat) {
  ggplot(node_stats, aes(x = .data[[stat]], fill = temperature, color = temperature)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    labs(x = stat, y = "Density",
         fill = "Temperature", color = "Temperature")
}
require(gridExtra)
norm_degree <- densityplot(stats_mb, "norm_degree")
norm_betweenness <- densityplot(stats_mb, "norm_betweenness")
norm_closeness <- densityplot(stats_mb, "norm_closeness")
eigenvector <- densityplot(stats_mb, "eigenvector")
grid.arrange(norm_degree, norm_betweenness, norm_closeness, eigenvector)


violinplot <- function(node_stats, stat) {
  ggplot(node_stats, aes(x = temperature, y = .data[[stat]], fill = temperature)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    theme_classic() +
    labs(x = "Temperature", y = stat, fill = "Temperature")
}
norm_degree <- violinplot(stats_mb, "norm_degree")
norm_betweenness <- violinplot(stats_mb, "norm_betweenness")
norm_closeness <- violinplot(stats_mb, "norm_closeness")
eigenvector <- violinplot(stats_mb, "eigenvector")
grid.arrange(norm_degree, norm_betweenness, norm_closeness, eigenvector)

 
