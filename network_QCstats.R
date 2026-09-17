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
  # build the path
  edges <- read.csv(file.path(dir, folder, "global", paste0("edges_", file_type, ".csv")))
  nodes <- read.csv(file.path(dir, folder, "global", paste0("nodes_", file_type, ".csv")))
  graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
}

networks <- list()

for (type in types) {
  for (qc in qcs) {     # make the graph for every network
    networks[[paste0(type, "_", qc)]] <- make_network(qc, type)
  }
}

# ------------------------------------------------------
# Node statistics distribution

# calculate degree, betweeennes, closeness and eigenvector
get_node_stats <- function(g, network_name) {
  tibble(network = network_name,
         network_type = sub("_.*", "", network_name),
         QC = sub(".*_", "", network_name),
         node = V(g)$name,
         degree = degree(g),
         norm_degree = degree(g, normalized = TRUE),
         betweenness = betweenness(g, weights = NA),
         norm_betweenness = betweenness(g, weights = NA, normalized = TRUE),
         closeness = closeness(g, weights = NA),
         norm_closeness = closeness(g, weights = NA, normalized = TRUE),
         eigenvector = eigen_centrality(g, weights = NA)$vector)
}    # all in the same table 

node_stats <- purrr::imap_dfr(
  networks,
  get_node_stats 
)
# now we need to change the type and QC to factor to graph
node_stats <- node_stats %>% mutate(network_type = factor(network_type, levels = c("MAG-MAG", "MAG-BGC", "MAG-MAG-rec")), 
                                    QC = factor(QC, levels = c("0", "08", "15")))
# BOXPLOT
boxplot <- function(node_stats, stat) {
  ggplot(node_stats, aes(x = QC, y = .data[[stat]], fill = network_type)) +
    geom_boxplot() +
    theme_minimal() +
    labs(x = "QC", y = stat, fill = "Network")
}
require(gridExtra)
norm_degree <- boxplot(node_stats, "norm_degree")
degree <- boxplot(node_stats, "degree")
betweenness <- boxplot(node_stats, "betweenness")
norm_betweenness <- boxplot(node_stats, "norm_betweenness")
closeness <- boxplot(node_stats, "closeness")
eigenvector <- boxplot(node_stats, "eigenvector")
grid.arrange(degree, norm_degree, norm_betweenness, betweenness, closeness, eigenvector)


# QQPLOT
degree_mm <- ggplot(node_stats %>% filter(network_type == "MAG-MAG", QC == "15"), aes(sample = norm_degree)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal() +
  labs(title = "MAG-MAG node degree")
degree_mb <- ggplot(node_stats %>% filter(network_type == "MAG-BGC", QC == "15"), aes(sample = norm_degree)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal() +
  labs(title = "MAG-BGC node degree")
grid.arrange(degree_mb, degree_mm, ncol = 2, nrow = 1)


# -----------------------------------------------------
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
  labs(title = "MAG-BGC networks", x = NULL, y = NULL, fill = "Jaccard")
grid.arrange(mag_bgc, mag_mag, mag_mag_r, nrow = 1, ncol = 3)


# ---------------------------------
# Shared edges between networks

library(VennDiagram)

# make edges unique identifiers 
make_edge_id <- function(data, node1, node2) { 
  apply( data[, c(node1, node2)], 1, function(x) { 
    paste(sort(x), collapse = "--") } ) 
}

# MAG-MAG
potentials_mm <- read.csv("2026-09-interactions/mOTUs_Species_Cluster/global/all_cases.csv") %>% 
  filter(oc_sites >= 1)    # de todas las combinaciones posibles, la que co-ocurren al menos una vez
lowq_mm <- read.csv("2026-09-interactions/mOTUs_Species_Cluster/global/oc_filt.csv") # estadisticamente significativas pero con sitios de baja calidad
highq_mm <- read.csv("2026-09-interactions08/mOTUs_Species_Cluster/global/oc_filt.csv") # con sitios de alta calidad
potentials_mm$edge <- make_edge_id(potentials_mm, "MAGi", "MAGj" )
lowq_mm$edge <- make_edge_id(lowq_mm, "MAGi", "MAGj" )
highq_mm$edge <- make_edge_id(highq_mm, "MAGi", "MAGj" )
# hacer un set
edges_mm <- list(Potential = unique(potentials_mm$edge), 
                 Low_QC = unique(lowq_mm$edge), 
                 High_QC = unique(highq_mm$edge))
myCol <- c("#56B4E9", "#E69F00", "#009E73")

# aristas que comparten
venn.diagram(
  x = list(Potential = edges_mm$Potential, Low_QC = edges_mm$Low_QC, High_QC = edges_mm$High_QC),
  category.names = c("Potential", "Low QC", "High QC"),
  filename = "MAG-MAG_venn.png",
  output = TRUE,
  # Output
  imagetype = "png", height = 480, width = 480, resolution = 300, compression = "lzw",
  # Circles
  lwd = 2, lty = "blank", fill = myCol,
  # Numbers
  cex = 0.6, fontfamily = "sans",
  # Set names
  cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans", rotation = 1
)


# MAG-BGC
potentials_mb <- read.csv("2026-09-interactions/mOTUs_Species_Cluster_gcc/global/all_cases.csv") %>%
  filter(oc_sites >= 1)
lowq_mb <- read.csv("2026-09-interactions/mOTUs_Species_Cluster_gcc/global/oc_filt.csv") # estadisticamente significativas pero con sitios de baja calidad
highq_mb <- read.csv("2026-09-interactions08/mOTUs_Species_Cluster_gcc/global/oc_filt.csv") # con sitios de alta calidad

potentials_mb$edge <- make_edge_id(potentials_mb, "Mags", "Bgcs" )
lowq_mb$edge <- make_edge_id(lowq_mb, "Mags", "Bgcs" )
highq_mb$edge <- make_edge_id(highq_mb, "Mags", "Bgcs" )

edges_mb <- list(Potential = unique(potentials_mb$edge), 
                 Low_QC = unique(lowq_mb$edge), 
                 High_QC = unique(highq_mb$edge))

venn.diagram(
  x = list(Potential = edges_mb$Potential, Low_QC = edges_mb$Low_QC, High_QC = edges_mb$High_QC),
  category.names = c("Potential", "Low QC", "High QC"),
  filename = "MAG-BGC_venn.png",
  output = TRUE,
  # Output
  imagetype = "png", height = 480, width = 480, resolution = 300, compression = "lzw",
  # Circles
  lwd = 2, lty = "blank", fill = myCol,
  # Numbers
  cex = 0.6, fontfamily = "sans",
  # Set names
  cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans", rotation = 1
)


# otra forma 
ggVennDiagram(
  edges_mb,
  label_alpha = 0,
  set_color = c("#56B4E9", "#E69F00", "#009E73")
) +
  scale_fill_gradient(
    low = "white",
    high = "white"
  ) +
  labs(
    title = "MAG-BGC interactions"
  )

# MAG-MAG-rec ???



