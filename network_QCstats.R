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

qcs <- c("0", "08", "15")
types <- c("MAG-MAG", "MAG-BGC", "MAG-MAG-rec")

make_network <- function(qc, type) {
  prefix <- ifelse(qc == "0", "", qc)
  dir <- paste0("2026-09-interactions", prefix)
  
  folder <- ifelse(type == "MAG-MAG", "mOTUs_Species_Cluster", "mOTUs_Species_Cluster_gcc")
  file_type <- ifelse(type == "MAG-BGC", "mb", "mm")
  
  edges <- read.csv(file.path(dir, folder, "global", paste0("edges_", file_type, ".csv")))
  nodes <- read.csv(file.path(dir, folder, "global", paste0("nodes_", file_type, ".csv")))
  graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
}

networks <- list()

for (type in types) {
  for (qc in qcs) {
    networks[[paste0(type, "_", qc)]] <- make_network(qc, type)
  }
}


# Node statistics distribution




# JACCARD
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

# graph heatmap
jaccard_to_df <- function(mat) {    # convestirlo a data frame
  as.data.frame(mat) %>%
    rownames_to_column("QC_1") %>%
    pivot_longer(cols = -QC_1,names_to = "QC_2", values_to = "Jaccard")
}
jaccard_mm_df <- jaccard_to_df(jaccard_mm)
jaccard_mb_df <- jaccard_to_df(jaccard_mb)
jaccard_mmr_df <- jaccard_to_df(jaccard_mmr)

mag_mag <- ggplot(jaccard_mm_df, aes(x = QC_2, y = QC_1, fill = Jaccard)) +
  geom_tile() +
  geom_text( aes(label = round(Jaccard, 2))) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  coord_equal() +
  theme_minimal() +
  labs(title = "MAG-MAG networks", x = NULL, y = NULL, fill = "Jaccard")
mag_mag_r <- ggplot(jaccard_mmr_df, aes(x = QC_2, y = QC_1, fill = Jaccard)) +
  geom_tile() +
  geom_text( aes(label = round(Jaccard, 2))) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  coord_equal() +
  theme_minimal() +
  labs(title = "MAG-MAG reconstructed networks", x = NULL, y = NULL, fill = "Jaccard")
mag_bgc <- ggplot(jaccard_mb_df, aes(x = QC_2, y = QC_1, fill = Jaccard)) +
  geom_tile() +
  geom_text( aes(label = round(Jaccard, 2))) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  coord_equal() +
  theme_minimal() +
  labs(title = "MAG-BGC reconstructed networks", x = NULL, y = NULL, fill = "Jaccard")
grid.arrange(mag_bgc, mag_mag, mag_mag_r, nrow = 1, ncol = 3)


# Shared edges between networks


