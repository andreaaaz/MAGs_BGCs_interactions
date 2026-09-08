################################################
### Network Statistics QC ######################
## Andrea Zermeño Díaz #########################
# september-2026 ###############################

# Compare the networks from di

#libraries
suppressPackageStartupMessages(library(tidyverse))
library(optparse)
library(purrr)
suppressPackageStartupMessages(library(igraph))

# arguments
# loading data
setwd("~/interactions/")

# mm: MAG-MAG, mb: MAG-BGC, mmr: MAG-MAG reconstructed, e: edges, n: nodes, NUMBERS: QC

mm_15_e <- read.csv("2026-09-interactions15/mOTUs_Species_Cluster/global/edges_mm.csv")
mm_15_n <- read.csv("2026-09-interactions15/mOTUs_Species_Cluster/global/nodes_mm.csv")
mb_15_e <- read.csv("2026-09-interactions15/mOTUs_Species_Cluster_gcc/global/edges_mb.csv")
mb_15_n <-  read.csv("2026-09-interactions15/mOTUs_Species_Cluster_gcc/global/nodes_mb.csv")
mmr_15_e <- read.csv("2026-09-interactions15/mOTUs_Species_Cluster_gcc/global/edges_mm.csv")
mmr_15_n <- read.csv("2026-09-interactions15/mOTUs_Species_Cluster_gcc/global/nodes_mm.csv")

mm_0_e <- read.csv("2026-09-interactions/mOTUs_Species_Cluster/global/edges_mm.csv")
mm_0_n <- read.csv("2026-09-interactions/mOTUs_Species_Cluster/global/nodes_mm.csv")
mb_0_e <- read.csv("2026-09-interactions/mOTUs_Species_Cluster_gcc/global/edges_mb.csv")
mb_0_n <-  read.csv("2026-09-interactions/mOTUs_Species_Cluster_gcc/global/nodes_mb.csv")
mmr_0_e <- read.csv("2026-09-interactions/mOTUs_Species_Cluster_gcc/global/edges_mm.csv")
mmr_0_n <- read.csv("2026-09-interactions/mOTUs_Species_Cluster_gcc/global/nodes_mm.csv")

mm_8_e <- read.csv("2026-09-interactions8/mOTUs_Species_Cluster/global/edges_mm.csv")
mm_8_n <- read.csv("2026-09-interactions8/mOTUs_Species_Cluster/global/nodes_mm.csv")
mb_8_e <- read.csv("2026-09-interactions8/mOTUs_Species_Cluster_gcc/global/edges_mb.csv")
mb_8_n <-  read.csv("2026-09-interactions8/mOTUs_Species_Cluster_gcc/global/nodes_mb.csv")
mmr_8_e <- read.csv("2026-09-interactions8/mOTUs_Species_Cluster_gcc/global/edges_mm.csv")
mmr_8_n <- read.csv("2026-09-interactions8/mOTUs_Species_Cluster_gcc/global/nodes_mm.csv")

# MAKE NETWORKS
make_network <- function(edges, nodes, directed = FALSE) {
  graph_from_data_frame(d = edges, vertices = nodes, directed = directed)
}
g_mm_0  <- make_network(mm_0_e, mm_0_n)
g_mm_8  <- make_network(mm_8_e, mm_8_n)
g_mm_15 <- make_network(mm_15_e, mm_15_n)

g_mb_0  <- make_network(mb_0_e, mb_0_n)
g_mb_8  <- make_network(mb_8_e, mb_8_n)
g_mb_15 <- make_network(mb_15_e, mb_15_n)

g_mmr_0  <- make_network(mmr_0_e, mmr_0_n)
g_mmr_8  <- make_network(mmr_8_e, mmr_8_n)
g_mmr_15 <- make_network(mmr_15_e, mmr_15_n)

# network list
networks <- list(
  "MAG-MAG_0"  = g_mm_0,
  "MAG-MAG_8"  = g_mm_8,
  "MAG-MAG_15" = g_mm_15,
  
  "MAG-BGC_0"  = g_mb_0,
  "MAG-BGC_8"  = g_mb_8,
  "MAG-BGC_15" = g_mb_15,
  
  "MAG-MAG-rec_0"  = g_mmr_0,
  "MAG-MAG-rec_8"  = g_mmr_8,
  "MAG-MAG-rec_15" = g_mmr_15
)


# Network Statistics
get_network_stats <- function(g) {
  comp <- components(g)
  tibble(n_nodes = vcount(g), 
         n_edges = ecount(g), 
         density = edge_density(g), 
         diameter = diameter(g, directed = FALSE, weights = NA), # no contar p-values como distancia
         n_components = comp$no, 
         giant_component = max(comp$csize), 
         giant_component_prop = max(comp$csize) / vcount(g),
         transitivity = transitivity(g, type = "global"))
}
network_stats <- imap_dfr(networks, ~ get_network_stats(.x) %>%
                            mutate(network = .y)) %>% 
  select(network, everything())
network_stats <- network_stats %>% # reacomodar
  separate(network, into = c("network_type", "QC"), sep = "_")


# GRAPH
network_stats$QC <- factor(network_stats$QC, levels = c(0, 8, 15))
require(gridExtra)
nodes <- ggplot(network_stats, aes(x = QC, y = n_nodes, color = network_type, group = network_type)) +
  geom_point(size = 2) +
  geom_line() +
  theme_minimal()
edges <- ggplot(network_stats, aes(x = QC, y = n_edges, color = network_type, group = network_type)) +
  geom_point(size = 2) +
  geom_line() +
  theme_minimal()
density <- ggplot(network_stats, aes(x = QC, y = density, color = network_type, group = network_type)) +
  geom_point(size = 2) +
  geom_line() +
  theme_minimal()
diameter <- ggplot(network_stats, aes(x = QC, y = diameter, color = network_type, group = network_type)) +
  geom_point(size = 2) +
  geom_line() +
  theme_minimal()
degree <- ggplot(network_stats, aes(x = QC, y = mean_degree, color = network_type, group = network_type)) +
  geom_point(size = 2) +
  geom_line() +
  theme_minimal()
component <- ggplot(network_stats, aes(x = QC, y = giant_component, color = network_type, group = network_type)) +
  geom_point(size = 2) +
  geom_line() +
  theme_minimal()
transitivity <- ggplot(network_stats, aes(x = QC, y = transitivity, color = network_type, group = network_type)) +
  geom_point(size = 2) +
  geom_line() +
  theme_minimal()
grid.arrange(nodes, edges, density, diameter, degree, component, transitivity)

# compare networks with Jaccard

# calcular el indice de jaccard de aristas entre dos redes
jaccard_edges <- function(g1, g2) {
  edges1 <- apply(as_edgelist(g1), 1, function(x) {
    paste(sort(x), collapse = "--")
  })
  edges2 <- apply(as_edgelist(g2), 1, function(x) {
    paste(sort(x), collapse = "--")
  })
  intersection <- length(intersect(edges1, edges2))
  union <- length(union(edges1, edges2))
  if (union == 0) {
    return(NA)
  }
  intersection / union
}
# todas las comparaciones entre redes
jaccard_matrix <- function(network_list) {
  n <- length(network_list)
  mat <- matrix(NA, nrow = n, ncol = n)
  rownames(mat) <- names(network_list)
  colnames(mat) <- names(network_list)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mat[i, j] <- jaccard_edges(
        network_list[[i]],
        network_list[[j]])
    }
  }
  mat
}

mm_networks <- networks[grepl("^MAG-MAG_[0-9]+$", names(networks))]
mb_networks <- networks[grepl("^MAG-BGC_", names(networks))]
mmr_networks <- networks[grepl("^MAG-MAG-rec_", names(networks))]

jaccard_mm <- jaccard_matrix(mm_networks)
jaccard_mb <- jaccard_matrix(mb_networks)
jaccard_mmr <- jaccard_matrix(mmr_networks)

# Node distribution 