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
library(bipartite)

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

temp_cols <- c("high" = "#B30000", "mid" = "#FC8D59", "low" = "#41B6C4")
temp_order <- c("high", "mid", "low", "global")

#------------------------------------------------------------------------
# NODE STATS

dir.create(file.path(opt$outdir, "nodes"), showWarnings = FALSE, recursive = TRUE)

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
  mutate(log_norm_betweenness = log1p(norm_betweenness)) %>%
  mutate(temperature = factor(temperature, levels = temp_order))

densityplot <- function(df, stat, exclude_global = FALSE){
  if(exclude_global){df <- df %>% filter(temperature != "global")}

  ggplot(df, aes(x = .data[[stat]], color = temperature, fill = temperature, linetype = temperature)) +
    geom_density(aes(fill = ifelse(temperature == "global", NA, 
                                   as.character(temperature))), alpha = 0.2) +
    scale_color_manual(values = c(temp_cols, "global" = "black"), 
                       breaks = c("high", "mid", "low", "global"), drop = TRUE) +
    scale_fill_manual(values = c(temp_cols, "global" = NA), 
                      breaks = c("high", "mid", "low"), drop = TRUE) +
    scale_linetype_manual(values = c("high" = "solid", "mid" = "solid", "low" = "solid","global" = "dashed"),
                          breaks = c("high", "mid", "low", "global"), drop = TRUE) +
    theme_classic() +
    labs(title = opt$type, x = stat, y = "Density", 
         color = "Temperature", fill = "Temperature",linetype = "Temperature") +
    guides(fill = "none", linetype = "none", color = guide_legend(override.aes = list(fill = NA, linetype = "solid")))
}

norm_degree <- densityplot(stats, "norm_degree")
norm_betweenness <- densityplot(stats, "log_norm_betweenness", exclude_global = TRUE)
norm_closeness <- densityplot(stats, "norm_closeness")
eigenvector <- densityplot(stats, "eigenvector")

# saving image and table
png(file.path(opt$outdir, paste0(opt$type, "nodes", "n_density.png")),width = 2000,height = 1600, res = 200) 
grid.arrange(norm_degree, norm_betweenness, norm_closeness, eigenvector)
dev.off()

write.csv(stats, file.path(opt$outdir, paste0(opt$type, "nodes", "n_stats.csv")), row.names = FALSE)

#---------------------------------------------------------------------
# EDGES STATS

dir.create(file.path(opt$outdir, "edges"), showWarnings = FALSE, recursive = TRUE)

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
edge_stats <- edge_stats %>%
  mutate(log_edge_betweenness = log1p(edge_betweenness)) # sale muy cargado al 0

edge_densityplot <- function(df, x, title){
  ggplot(df, aes(x = .data[[x]], color = temperature, linetype = temperature)) +
    geom_density(aes(fill = ifelse(temperature == "global", NA, as.character(temperature))),alpha = 0.2) +
    scale_color_manual(values = c(temp_cols, "global" = "black"), 
                       breaks = c("high", "mid", "low", "global"), drop = TRUE) +
    scale_fill_manual(values = c(temp_cols, "global" = NA), 
                      breaks = c("high", "mid", "low"),drop = TRUE) +
    scale_linetype_manual(values = c("high" = "solid", "mid" = "solid", "low" = "solid", "global" = "dashed"),
                          breaks = c("high", "mid", "low", "global"), drop = TRUE) +
    theme_classic() +
    labs(title = title, x = x, y = "Density", color = "Temperature", fill = "Temperature", linetype = "Temperature") +
    guides(fill = "none", linetype = "none", color = guide_legend(override.aes = list(fill = NA, linetype = "solid")))
}

edge_betweenness_plot <- edge_densityplot(edge_stats, "log_edge_betweenness", paste(opt$type, "- Edge betweenness")) +
  labs(x = "log(1 + edge betweenness)")
edge_connectivity_plot <- edge_densityplot(edge_stats, "edge_connectivity", paste(opt$type, "- Edge connectivity"))

png(file.path(opt$outdir, paste0(opt$type, "e_density.png")), width = 2000, height = 800, res = 200)
grid.arrange( edge_betweenness_plot,edge_connectivity_plot,ncol = 2)
dev.off()

write.csv(edge_stats, file.path(opt$outdir, paste0(opt$type, "e_stats.csv")), row.names = FALSE) # necesito que se guarde en temp_stats/edges/

#---------------------------------------------------------------------
# DISSIMILARITY OF SPECIES

dir.create(file.path(opt$outdir, "dissimilarity"), showWarnings = FALSE, recursive = TRUE)

# calc dissimilarity
edge_keys <- function(g) {
  el <- as_edgelist(g, names = TRUE)
  data.frame(a = pmin(el[,1], el[,2]), b = pmax(el[,1], el[,2]), stringsAsFactors = FALSE)
}

sorensen <- function(x, y) {           # x, y = vectores de aristas
  a <- length(intersect(x, y))
  b <- length(setdiff(x, y))
  c <- length(setdiff(y, x))
  if ((2*a + b + c) == 0) return(NA_real_)
  (b + c) / (2*a + b + c)
}

beta_links <- function(g1, g2) {
  shared <- intersect(V(g1)$name, V(g2)$name)
  
  k1 <- edge_keys(g1); k2 <- edge_keys(g2)
  key <- function(k) paste(k$a, k$b, sep = "||")
  only_shared <- function(k) k[k$a %in% shared & k$b %in% shared, ]
  
  WN <- sorensen(key(k1), key(k2))                              # total
  OS <- sorensen(key(only_shared(k1)), key(only_shared(k2)))    # rewiring
  S  <- sorensen(V(g1)$name, V(g2)$name)                        # especies
  data.frame(S = S, 
             OS = OS, 
             ST = WN - OS, 
             WN = WN, 
             prop_rewiring = OS / WN)
}

# todas las parejas de temperaturas
pairs <- combn(names(networks), 2, simplify = FALSE)
beta_df <- do.call(rbind, lapply(pairs, function(p) {
  cbind(t1 = p[1], t2 = p[2], beta_links(networks[[p[1]]], networks[[p[2]]]))
}))
write.csv(beta_df,file.path(opt$outdir, paste0(opt$type, "dissimilarity" ,"d_stats.csv")), row.names = FALSE)

# graph 
beta_vs <- beta_df %>%
  mutate(comparison = paste(t1, t2, sep = " vs "),
         comparison = factor(comparison, levels = c("global vs low", "global vs mid", "global vs high",
                                               "low vs mid", "low vs high", "mid vs high")))
beta_long <- beta_vs %>%
  select(comparison, OS, ST) %>%
  pivot_longer(cols = c(OS, ST), names_to = "component", values_to = "value") %>%
  mutate(component = recode(component, OS = "Rewiring", ST = "Turnover"))

dissimilarity <- ggplot(beta_long, aes(x = comparison, y = value, fill = component)) +
  geom_col(width = 0.7) +
  geom_text(data = beta_vs, aes(x = comparison, y = WN, label = round(WN, 2)),
    inherit.aes = FALSE, vjust = -0.5, size = 3.5) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
  labs(x = "Temperature comparison", y = "Network dissimilarity (βWN)",
    fill = "Component", title = "Decomposition of network dissimilarity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

ggsave(dissimilarity, )

#---------------------------------------------------------------------
# NESTEDNESS AND MODULARITY

dir.create(file.path(opt$outdir, "nest_mod"), showWarnings = FALSE, recursive = TRUE)

# matriz de incidencia 
to_web <- function(g, type_col = "type") {
  tp <- vertex_attr(g, type_col)
  V(g)$type <- tp == unique(tp)[1]
  stopifnot(is_bipartite(g))
  
  web <- as_biadjacency_matrix(g, types = V(g)$type, sparse = FALSE)
  web[rowSums(web) > 0, colSums(web) > 0, drop = FALSE]
}
net_metrics <- function(g) {
  web <- to_web(g)
  # Nestedness
  nodf <- unname(networklevel(web, index = "NODF"))
  # Modularity
  mod <- computeModules(web)
  Q <- mod@likelihood
  data.frame(n_rows = nrow(web),
             n_cols = ncol(web), 
             n_links = sum(web > 0),
             NODF = nodf, Q = Q)
}

metrics <- do.call(rbind, lapply(names(networks), function(t) {
  message("Procesando: ", t)
  cbind(temp = t, net_metrics(networks[[t]]))
}))

write.csv()




# load data out of nextflow

# base_dir <- "~/interactions/08_interactions/"
# base_dirs <- c(
#   mm   = "mOTUs_Species_Cluster",
#   mmr  = "mOTUs_Species_Cluster_gcc",
#   mb   = "mOTUs_Species_Cluster_gcc",
#   cons = "mOTUs_Species_Cluster_gcc"
# )
# # Network type to analyze
# opt$type <- "mb"  # mm, mmr, mb, cons
# temps <- c("global", "low", "mid", "high")
# # Function to create networks
# make_network <- function(edges_file, nodes_file) {
#   edges <- read.csv(edges_file, stringsAsFactors = FALSE)
#   nodes <- read.csv(nodes_file, stringsAsFactors = FALSE)
#   if ("node" %in% names(nodes))
#     nodes <- nodes %>% rename(id = node)
#   graph_from_data_frame(d = edges, vertices = nodes,directed = FALSE)
# }
# 
# # Load networks
# networks <- purrr::map(temps, function(temp) {
#   current_dir <- file.path(base_dir, base_dirs[opt$type], temp)
#   edges_file <- file.path(current_dir, paste0("edges_", opt$type, ".csv"))
#   nodes_file <- file.path(current_dir, paste0("nodes_", opt$type, ".csv"))
#   make_network(edges_file, nodes_file)
# }) %>%
#   set_names(temps)
# 

