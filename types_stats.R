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
