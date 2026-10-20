################################################
### Network stats by type  ######################
## Andrea Zermeño Díaz #########################
# october-2026 ###############################

# compare the type of networks

# libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(tidyverse))
library(ggplot2)

# file args 
option_list <- list(make_option(c("-t", "--temp"), type="character"),
                    make_option(c("--edges_mb"), type="character"),
                    make_option(c("--nodes_mb"), type="character"),
                    make_option(c("--edges_mmr"), type="character"),
                    make_option(c("--nodes_mmr"), type="character"),
                    make_option(c("--edges_mm"), type="character"),
                    make_option(c("--nodes_mm"), type="character"),
                    make_option(c("--edges_cons"), type="character"),
                    make_option(c("--nodes_cons"), type="character"),
                    make_option(c("-o", "--outdir"), type="character"))
opt <- parse_args(OptionParser(option_list = option_list))

# network load
make_network <- function(edges_file, nodes_file, directed = FALSE){
  edges <- read.csv(edges_file)
  nodes <- read.csv(nodes_file)
  
  if("node" %in% names(nodes)) # para consensus
    nodes <- nodes %>% rename(id = node)
  
  graph_from_data_frame(d = edges, vertices = nodes, directed = directed)
}

networks <- list("MAG-BGC" = make_network(opt$edges_mb,   opt$nodes_mb),
                 "MAG-MAG-rec" = make_network(opt$edges_mmr,  opt$nodes_mmr),
                 "MAG-MAG" = make_network(opt$edges_mm,   opt$nodes_mm),
                 "consensus" = make_network(opt$edges_cons, opt$nodes_cons))

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

node_stats <- purrr::imap_dfr(networks, ~get_node_stats(.x, .y, opt$temp))

densityplot <- function(node_stats, stat) {
  ggplot(node_stats, aes(x = .data[[stat]], fill = network_type, color = network_type)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    labs(x = stat, y = "Density",
         fill = "Network type", color = "Network type")
}

norm_degree <- densityplot(stats_global, "norm_degree") 
norm_betweenness <- densityplot(stats_global, "norm_betweenness")
norm_closeness <- densityplot(stats_global, "norm_closeness")
eigenvector <- densityplot(stats_global, "eigenvector")

png(opt$outdir, width = 2000, height = 1600, res = 200)
grid.arrange(norm_degree, norm_betweenness, norm_closeness, eigenvector, ncol = 2)
dev.off()

