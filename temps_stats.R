################################################
### Network stats by temperature  ######################
## Andrea Zermeño Díaz #########################
# october-2026 ###############################

# compare the effect of the temperature on the networks

# libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gridExtra))
library(ggplot2)

# file args 
option_list <- list(make_option("--type", type = "character", help = "Network type"),
                    make_option("--temps", type = "character", help = "Comma-separated temperatures"),
                    make_option("--edges", type = "character", help = "Space-separated paths to edge files"),
                    make_option("--nodes", type = "character", help = "Space-separated paths to node files"),
                    make_option("--outdir", type = "character", default = ".", help = "Output directory"))
opt <- parse_args(OptionParser(option_list = option_list))

# parse input list
temps <- strsplit(opt$temps, ",", fixed = TRUE)[[1]]
edges_files <- strsplit(opt$edges, ",", fixed = TRUE)[[1]]
nodes_files <- strsplit(opt$nodes, ",", fixed = TRUE)[[1]]

edges_files <- trimws(edges_files)
nodes_files <- trimws(nodes_files)

# make networks
make_network <- function(edges_file, nodes_file) {
  edges <- read.csv(edges_file, stringsAsFactors = FALSE)
  nodes <- read.csv(nodes_file, stringsAsFactors = FALSE)
  
  if("node" %in% names(nodes)) # para consensus
    nodes <- nodes %>% rename(id = node)
  
  graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
}

networks <- purrr::map2(edges_files, nodes_files, make_network)
names(networks) <- temps

#------------------------------------------------------------------------
# NODE STATS
# calculate degree, betweeennes, closeness and eigenvector
get_node_stats <- function(g, network_type, temp) {
  tibble(network_type = network_type,
         temperature = temp,
         node = V(g)$name,
         degree = degree(g),
         norm_degree = degree(g, normalized = TRUE),
         betweenness = betweenness(g, weights = NA),
         norm_betweenness = betweenness(g, weights = NA, normalized = TRUE),
         closeness = closeness(g, weights = NA),
         norm_closeness = closeness(g, weights = NA, normalized = TRUE),
         eigenvector = eigen_centrality(g, weights = NA)$vector)
}

stats <- purrr::imap_dfr(networks,~ get_node_stats(.x, opt$type, .y))
stats <- stats %>%
  mutate(log_norm_betweenness = log1p(norm_betweenness))

densityplot <- function(df, stat){
  ggplot(df, aes(x = .data[[stat]], fill = temperature, color = temperature)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    labs(title = opt$type, x = stat, y = "Density", color = "Temperature", fill = "Temperature")
}

norm_degree <- densityplot(stats, "norm_degree")
norm_betweenness <- densityplot(stats, "log_norm_betweenness")
norm_closeness <- densityplot(stats, "norm_closeness")
eigenvector <- densityplot(stats, "eigenvector")

# saving image and table
png(file.path(opt$outdir, paste0(opt$type, ".png")),width = 2000,height = 1600, res = 200)
grid.arrange(norm_degree, norm_betweenness, norm_closeness, eigenvector)
dev.off()

write.csv(stats, file.path(opt$outdir, paste0(opt$type, "_stats.csv")), row.names = FALSE)

#---------------------------------------------------------------------
# EDGES STATS
get_edge_stats <- function(g, network_type, temp) {
  
  edge_ends <- ends(g, E(g), names = TRUE) # Identificadores de los extremos de cada arista
  eb <- edge_betweenness(g, directed = FALSE, weights = NA) # Centralidad de intermediación de cada arista
  
  # Conectividad local entre los extremos de cada arista
  ec <- vapply(seq_len(ecount(g)), function(i) {
    ends_i <- edge_ends[i, ]
    edge_connectivity(g, source = ends_i[1], target = ends_i[2])
  }, numeric(1))
  
  tibble(network_type = network_type,
         temperature = temp,
         source = edge_ends[, 1],
         target = edge_ends[, 2],
         edge_betweenness = as.numeric(eb),
         edge_connectivity = as.numeric(ec))
}

edge_stats <- purrr::imap_dfr(networks, ~ get_edge_stats(.x, opt$type, .y))


edge_betweenness_plot <- ggplot(edge_stats, aes(x = log1p(edge_betweenness), color = temperature, fill = temperature)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  labs(title = paste(opt$type, "- Edge betweenness"), x = "log(1 + edge betweenness)", y = "Density")
edge_connectivity_plot <- ggplot(edge_stats,aes(x = edge_connectivity, color = temperature, fill = temperature)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  labs( title = paste(opt$type, "- Edge connectivity"), x = "Edge connectivity", y = "Density")

png(file.path(opt$outdir, paste0(opt$type, "_edge_stats.png")), width = 2000, height = 800, res = 200)
grid.arrange( edge_betweenness_plot,edge_connectivity_plot,ncol = 2)
dev.off()

write.csv(edge_stats, file.path(opt$outdir, paste0(opt$type, "_edge_stats.csv")), row.names = FALSE)

# DISSIMILARITY OF SPECIES
