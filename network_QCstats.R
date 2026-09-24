################################################
### Network Statistics QC ######################
## Andrea Zermeño Díaz #########################
# september-2026 ###############################

# Compare the networks from di

#libraries
suppressPackageStartupMessages(library(tidyverse))
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
mag_mag <- node_stats %>%
  filter(network_type == "MAG-MAG", QC %in% c("0", "08", "15"))
mag_bgc <- node_stats %>%
  filter(network_type == "MAG-BGC", QC %in% c("0", "08", "15"))
mag_mag_rec <- node_stats %>%
  filter(network_type == "MAG-MAG-rec", QC %in% c("0", "08", "15"))

q <- seq(0, 1, length.out = 100)

qq_mm <- tibble(
  QC_0 = quantile(mag_mag$norm_degree[mag_mag$QC == "0"], probs = q, na.rm = TRUE),
  QC_08 = quantile(mag_mag$norm_degree[mag_mag$QC == "08"], probs = q, na.rm = TRUE),
  QC_15 = quantile(mag_mag$norm_degree[mag_mag$QC == "15"], probs = q, na.rm = TRUE)
)
qq_mb <- tibble(
  QC_0 = quantile(mag_bgc$norm_degree[mag_bgc$QC == "0"], probs = q, na.rm = TRUE),
  QC_08 = quantile(mag_bgc$norm_degree[mag_bgc$QC == "08"], probs = q, na.rm = TRUE),
  QC_15 = quantile(mag_bgc$norm_degree[mag_bgc$QC == "15"], probs = q, na.rm = TRUE)
)
qq_mmr <- tibble(
  QC_0 = quantile(mag_mag_rec$norm_degree[mag_mag_rec$QC == "0"], probs = q, na.rm = TRUE),
  QC_08 = quantile(mag_mag_rec$norm_degree[mag_mag_rec$QC == "08"], probs = q, na.rm = TRUE),
  QC_15 = quantile(mag_mag_rec$norm_degree[mag_mag_rec$QC == "15"], probs = q, na.rm = TRUE)
)

plot_mm <- ggplot(qq_mm, aes(x = QC_0)) + 
  geom_point(aes(y = QC_08, color = "QC 08")) +
  geom_point(aes(y = QC_15, color = "QC 15")) +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  labs(x = "Quantiles QC 0", y = "Quantiles", color = "QC", title = "MAG-MAG normalized degree")
plot_mb <- ggplot(qq_mb, aes(x = QC_0)) + 
  geom_point(aes(y = QC_08, color = "QC 08")) +
  geom_point(aes(y = QC_15, color = "QC 15")) +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  labs(x = "Quantiles QC 0", y = "Quantiles", color = "QC", title = "MAG-BGC normalized degree")
plot_mmr <- ggplot(qq_mmr, aes(x = QC_0)) + 
  geom_point(aes(y = QC_08, color = "QC 08")) +
  geom_point(aes(y = QC_15, color = "QC 15")) +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  labs(x = "Quantiles QC 0", y = "Quantiles", color = "QC", title = "MAG-MAG-rec normalized degree")
grid.arrange(plot_mm, plot_mb, plot_mmr, ncol = 3, nrow = 1)


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
library(eulerr)
library(UpSetR)
library(ggVennDiagram)
library(VennDiagram)

# make edges unique identifiers 
make_edge_id <- function(data, node1, node2) { 
  apply( data[, c(node1, node2)], 1, function(x) { 
    paste(sort(x), collapse = "--") } ) 
}


# MAG-MAG -----------
potentials_mm <- read.csv("potentials_mm.csv") 
lowq_mm <- read.csv("2026-09-interactions/mOTUs_Species_Cluster/global/oc_filt.csv") # estadisticamente significativas pero con sitios de baja calidad
highq_mm <- read.csv("2026-09-interactions08/mOTUs_Species_Cluster/global/oc_filt.csv") # con sitios de alta calidad
lowq_mmr <- read.csv("2026-09-interactions/mOTUs_Species_Cluster_gcc/global/edges_mm.csv") # mam-mag recontruida
highq_mmr <- read.csv("2026-09-interactions08/mOTUs_Species_Cluster_gcc/global/edges_mm.csv")


potentials_mm$edge <- make_edge_id(potentials_mm, "MAGi", "MAGj" )
lowq_mm$edge <- make_edge_id(lowq_mm, "MAGi", "MAGj" )
highq_mm$edge <- make_edge_id(highq_mm, "MAGi", "MAGj" )
lowq_mmr$edge <- make_edge_id(lowq_mmr, "source", "target" )
highq_mmr$edge <- make_edge_id(highq_mmr, "source", "target" )

# hacer un UpSet
edges_global <- list(Potentials = unique(potentials_mm$edge),
                 Rec_Low_QC = unique(lowq_mmr$edge), 
                 Rec_High_QC = unique(highq_mmr$edge),
                 Low_QC = unique(lowq_mm$edge),
                 High_QC = unique(highq_mm$edge))
edges_global2 <- list(Rec_Low_QC = unique(lowq_mmr$edge), 
                 Rec_High_QC = unique(highq_mmr$edge),
                 Low_QC = unique(lowq_mm$edge),
                 High_QC = unique(highq_mm$edge))

upset(fromList(edges_global), 
      sets = c("Potentials", "Low_QC", "High_QC", "Rec_Low_QC","Rec_High_QC"), 
      keep.order = FALSE, order.by = "freq", mainbar.y.label = "Number of MAG-MAG edges",
      sets.x.label = "Number of edges")


# euler Diagram global
eulerr_global <- euler(list(Potentials = unique(potentials_mm$edge),
                         RecLow_QC = unique(lowq_mmr$edge), 
                         RecHigh_QC = unique(highq_mmr$edge),
                         Low_QC = unique(lowq_mm$edge),
                         High_QC = unique(highq_mm$edge)))
plot(eulerr_global, fills = TRUE, alpha = 0.8, quantities = TRUE, 
     labels = TRUE, edges = TRUE, proportional = FALSE, main = "Global MAG-MAG interactions")
# con zoom (ignora potenciales)
eulerr_global2 <- euler(list(RecLow_QC = unique(lowq_mmr$edge), 
                         RecHigh_QC = unique(highq_mmr$edge),
                         Low_QC = unique(lowq_mm$edge),
                         High_QC = unique(highq_mm$edge)))
plot(eulerr_global2, fills = TRUE, alpha = 0.8, quantities = TRUE,  
     labels = TRUE, edges = TRUE, proportional = TRUE, main = "Global MAG-MAG interactions")
# MAG-MAG
edges_mm <- list(Potentials = unique(potentials_mm$edge),
                  Low_QC = unique(lowq_mm$edge),
                  High_QC = unique(highq_mm$edge))
myCol <- c("#56B4E9", "#E69F00", "#009E73")
venn.diagram(x = edges_mm, category.names = c("Potentials", "Low QC", "High QC"), 
             filename = "MAG-MAG_venn.png", output = TRUE, # Output
             imagetype = "png", height = 480, width = 480, resolution = 300, compression = "lzw",
             # Circles
             lwd = 2, lty = "blank", fill = myCol,
             # Numbers
             cex = 0.6, fontfamily = "sans",
             # Set names
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans", rotation = 1)
# Reconstructed MAG-MAG
edges_mmr <- list(Potentials = unique(potentials_mm$edge),
                   RecLow_QC = unique(lowq_mmr$edge), 
                   RecHigh_QC = unique(highq_mmr$edge))
venn.diagram(x = edges_mmr, category.names = c("Potentials", "Low QC", "High QC"), 
             filename = "MAG-MAG-rec_venn.png", output = TRUE, # Output
             imagetype = "png", height = 480, width = 480, resolution = 300, compression = "lzw",
             # Circles
             lwd = 2, lty = "blank", fill = myCol,
             # Numbers
             cex = 0.6, fontfamily = "sans",
             # Set names
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans", rotation = 1)

# Other venn diagram
venn <- Venn(edges_global2)
data <- process_data(venn)
set_colors <- c("#F7AA14", "#F5D000", "#009E73", "#50C058")
ggplot() +3
  geom_polygon(
    aes(X, Y, group = id, fill = id, color = id), data = venn_setedge(data), alpha = 0.55, linewidth = 1) +
  scale_fill_manual(values = set_colors) +
  scale_color_manual(values = set_colors) +
  geom_text(aes(X, Y, label = name), data = venn_setlabel(data)) +
  # números de conteo por región
  geom_label(aes(X, Y, label = count), data = venn_regionlabel(data), fill = "white", 
             alpha = 0, label.size = NA, size = 3.5) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "MAG-MAG interactions")

# MAG-BGC ------------
potentials_mb <- read.csv("potentials_mb.csv") 
lowq_mb <- read.csv("2026-09-interactions/mOTUs_Species_Cluster_gcc/global/oc_filt.csv") # estadisticamente significativas pero con sitios de baja calidad
highq_mb <- read.csv("2026-09-interactions08/mOTUs_Species_Cluster_gcc/global/oc_filt.csv") # con sitios de alta calidad

potentials_mb$edge <- make_edge_id(potentials_mb, "Mags", "Bgcs" )
lowq_mb$edge <- make_edge_id(lowq_mb, "Mags", "Bgcs" )
highq_mb$edge <- make_edge_id(highq_mb, "Mags", "Bgcs" )

edges_mb <- list(Potential = unique(potentials_mb$edge), 
                 Low_QC = unique(lowq_mb$edge), 
                 High_QC = unique(highq_mb$edge))
# UpSet
upset(fromList(edges_mb), 
      sets = c("Potential", "Low_QC", "High_QC"), 
      keep.order = FALSE, order.by = "freq", mainbar.y.label = "Number of MAG-BGC edges",
      sets.x.label = "Number of edges")
# euler Diagram 
edges_mb <- list(Potentials = unique(potentials_mb$edge),
                  Low_QC = unique(lowq_mb$edge), 
                  High_QC = unique(highq_mb$edge))
venn.diagram(x = edges_mb, category.names = c("Potentials", "Low QC", "High QC"), 
             filename = "MAG-BGC_venn.png", output = TRUE, # Output
             imagetype = "png", height = 480, width = 480, resolution = 300, compression = "lzw",
             # Circles
             lwd = 2, lty = "blank", fill = myCol,
             # Numbers
             cex = 0.6, fontfamily = "sans",
             # Set names
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans", rotation = 1)

# MAG-MAG CONSENSUS

# este subset tiene las aristas que comparten las redes MAG-MAG y MAG-MAG reconstruida 
# en QC = 0 y QC = 8

highq_mmr <- highq_mmr %>%
  rename(MAGi = source, MAGj = target)
lowq_mmr <- lowq_mmr %>%
  rename(MAGi = source, MAGj = target)

consensus <- Reduce(
  function(x, y) inner_join(x, y, by = c("MAGi", "MAGj")),
  list(highq_mm, lowq_mm, highq_mmr, lowq_mmr)
)





