################################################
### Network Statistics QC ######################
## Andrea Zermeño Díaz #########################
# september-2026 ###############################

# Compare the networks from di

#libraries
suppressPackageStartupMessages(library(tidyverse))
library(optparse)
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
vcount(g)              # nodos
ecount(g)              # aristas
edge_density(g)        # densidad
mean(degree(g))        # grado promedio
components(g)$no       # componentes
diameter(g)            # diámetro
mean_distance(g)       # distancia promedio
transitivity(g)        # clustering/transitividad




