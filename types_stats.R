################################################
### Network stats by type  ######################
## Andrea Zermeño Díaz #########################
# october-2026 ###############################

# compare the type of networks

# libraries
library(igraph)
library(tidyverse)
library(ggplot2)

# data load
setwd("~/interactions/2026-09-interactions08/")
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

networks_global <- list()
networks_high <- list()
networks_mid <- list()
networks_low <- list()

for (type in types) {
  networks_global[[type]] <- make_network("global", type)
  networks_high[[type]] <- make_network("high", type)
  networks_mid[[type]] <- make_network("mid", type)
  networks_low[[type]] <- make_network("low", type)
}

# NODE STATS# calculate degree, betweeennes, closeness and eigenvector
get_node_stats <- function(g, network_name, network_type) {
  tibble(network = network_name,
         temperature = network_type,
         network_type = network_name,
         node = V(g)$name,
         degree = degree(g),
         norm_degree = degree(g, normalized = TRUE),
         betweenness = betweenness(g, weights = NA),
         norm_betweenness = betweenness(g, weights = NA, normalized = TRUE),
         closeness = closeness(g, weights = NA),
         norm_closeness = closeness(g, weights = NA, normalized = TRUE),
         eigenvector = eigen_centrality(g, weights = NA)$vector)
}

# global
stats_global <- purrr::imap_dfr(networks_global, ~ get_node_stats(.x, .y, "global"))
stats_global <- purrr::imap_dfr(networks_low, ~ get_node_stats(.x, .y, "low"))

densityplot <- function(node_stats, stat) {
  ggplot(node_stats, aes(x = .data[[stat]], fill = network_type, color = network_type)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    labs(x = stat, y = "Density",
         fill = "Network type", color = "Network type")
}

require(gridExtra)
norm_degree <- densityplot(stats_global, "norm_degree") 
norm_betweenness <- densityplot(stats_global, "norm_betweenness")
norm_closeness <- densityplot(stats_global, "norm_closeness")
eigenvector <- densityplot(stats_global, "eigenvector")
grid.arrange(norm_degree, norm_betweenness, norm_closeness, eigenvector)


# consensus
consensus_edges <- c(global = "~/interactions/consensus_edges.csv", 
                     high = "~/interactions/consensus_edges_h.csv",
                     mid = "~/interactions/consensus_edges_m.csv",
                     low = "~/interactions/consensus_edges_l.csv")

consensus_nodes <- c(global = "~/interactions/consensus_nodes.csv",
                     high = "~/interactions/consensus_nodes_h.csv",
                     mid = "~/interactions/consensus_nodes_m.csv",
                     low = "~/interactions/consensus_nodes_l.csv")
for (temperature in names(consensus_edges)) {
  edges <- read.csv(consensus_edges[temperature])
  nodes <- read.csv(consensus_nodes[temperature])
  g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  
  if (temperature == "global") {
    networks_global[["consensus"]] <- g
  } else if (temperature == "high") {
    networks_high[["consensus"]] <- g
  } else if (temperature == "mid") {
    networks_mid[["consensus"]] <- g
  } else if (temperature == "low") {
    networks_low[["consensus"]] <- g
  }
}





