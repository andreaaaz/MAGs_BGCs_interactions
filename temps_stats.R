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
suppressPackageStartupMessages(library(bipartite))
library(ggplot2)

# file args 
option_list <- list(make_option("--type", type = "character", help = "Network type"),
                    make_option("--temps", type = "character", help = "Comma-separated temperatures"),
                    make_option("--edges", type = "character", help = "Space-separated paths to edge files"),
                    make_option("--nodes", type = "character", help = "Space-separated paths to node files"),
                    make_option("--outfile", type = "character", help = "Output PNG file"),
                    make_option("--stats_out", type = "character", help = "Output CSV file"))
opt <- parse_args(OptionParser(option_list = option_list))

# parse input list

temps <- strsplit(opt$temps, ",", fixed = TRUE)[[1]]
edges_files <- strsplit(trimws(opt$edges), "\\s+")[[1]]
nodes_files <- strsplit(trimws(opt$nodes), "\\s+")[[1]]
# quitar cuando corra!!!!
if (length(temps) != length(edges_files) || length(temps) != length(nodes_files)) {
  stop("The number of temperatures, edge files and node files must match.")
}

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
png(opt$outfile,width = 2000,height = 1600, res = 200)
grid.arrange(norm_degree, norm_betweenness, norm_closeness, eigenvector)
dev.off()

write.csv(stats, file = opt$stats_out, row.names = FALSE)

# DISSIMILARITY OF SPECIES
