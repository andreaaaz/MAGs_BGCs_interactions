#######################################
#### NETWORKS MAGs-BGCs ###############
## Andrea Zermeño Díaz ################
# may-2026 ##########################

# libraries
suppressPackageStartupMessages(library(tidyverse))
library(ggplot2)
library(optparse)
suppressPackageStartupMessages(library(igraph))
library(paletteer)

### ARGS ---- 
option_list <- list(
  make_option(c("-m", "--microbial_lineage"), type="character", default="mOTUs_Species_Cluster", help="Name of the microbial lienage"),
  make_option(c("-b", "--bgc_groups"), type="character", default="gcc", help="Name of the grou"),
  make_option(c("-f", "--file_mb"), type="character", help="file with oc cases"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-i", "--indir"), type="character", help="Input directory"),
  make_option(c("-w", "--workdir"), type="character", help="Working directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
mag_lineage <- opt$microbial_lineage
bgc_group <- opt$bgc_groups

### DATA LOAD ----
# metadata
meta_mags <- read.csv(file = paste0(opt$indir, 'metadata.csv'), header = TRUE)
meta_bgcs <- read.csv(file = paste0(opt$indir, 'bgcs_metadata.csv'), header = TRUE)
# results
cases <- read.csv(opt$file_mb)
# functions
source(paste0(opt$workdir, "functions.R"))
# prep tables for the MAG-MAG filter
mags_by_sites <- prep_mags(meta_mags, mag_lineage)
bgcs_by_sites<- prep_bgcs(meta_bgcs, bgc_group)
# deberia de filtrar tambien por temperatura????

#-------------------------
#### MAG<-BGC NETWORK ####
#-------------------------

### NODES AND EDGES ----
nodes <- tibble(
  id = unique(c(cases$Mags, cases$Bgcs)),    
  type = ifelse(id %in% cases$Mags, "MAG", "BGC"))

edges <- cases %>%     #cases are already the edges, just rename it
  rename(source = Bgcs, target = Mags, weight = fdr_pval_o) 

### GRAPH AND DEGREE ----
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
nodes$degree <- degree(g)

### Add representative groups ----
rep_bgcs <- meta_bgcs %>%    # table with all the GCCs and their products
  group_by(.data[[bgc_group]], products) %>%   
  summarise(n = n(), .groups = "drop_last") %>%   # count the most common
  slice_max(order_by = n, n = 1) %>%             # between each GCC
  select(id = all_of(bgc_group), rep_bgc = products)
  
rep_mags <- meta_mags %>%   #table with all the mOTUs and their families
  group_by(.data[[mag_lineage]], family) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  slice_max(order_by = n, n = 1) %>%
  select(id = all_of(mag_lineage), rep_mag = family)
# add the representative groups to the node information
nodes <- nodes %>%
  left_join(rep_mags, by = "id") %>%
  left_join(rep_bgcs, by = "id") %>%
  mutate(color_group = if_else(type == "MAG", rep_mag, rep_bgc)) %>%
  select(-rep_bgc, -rep_mag)

# ### Top nodes ----
# # BUG: change to detect the top nodes without repetiton
# top_mag_nodes <- nodes %>%
#   filter(type == "MAG") %>%
#   arrange(desc(degree)) %>%
#   slice_head(n = 12)  # adapt the number to the colors in the palette
# top_bgc_nodes <- nodes %>%
#   filter(type == "BGC") %>% 
#   arrange(desc(degree)) %>%
#   slice_head(n = 12)
# top_mag_groups <- unique(top_mag_nodes$color_group)
# top_bgc_groups <- unique(top_bgc_nodes$color_group)
# 
# ## Colors for top nodes
# mag_colors <- substr(
#   as.vector(paletteer::paletteer_d(
#     palette = "ggthemes::calc",
#     n = length(top_mag_groups))), 1, 7) # adapt the number to the nodes
# bgc_colors <- substr(
#   as.vector(paletteer::paletteer_d(
#     palette = "ggthemes::gdoc",
#     n = length(top_bgc_groups))), 1, 7)
# 
# nodes <- nodes %>%
#   mutate(color = case_when(
#     type == "MAG" & color_group %in% top_mag_groups ~ 
#       mag_colors[match(color_group, top_mag_groups)],
#     type == "BGC" & color_group %in% top_bgc_groups ~ 
#       bgc_colors[match(color_group, top_bgc_groups)],
#     TRUE ~ "darkgrey"))


### Crear red en Cytoscape ----
# cytoscapePing()
# createNetworkFromDataFrames(
#  nodes = nodes,
#  edges = edges,
#  title = "oc mOTUs-GCCs",
#  collection = "Interacciones MAGs-BGCs")

### Sub y sobre abundancia productos BGCs  ----

# # esperados (dataset global)
# esp <- meta_bgcs %>%
#   count(products) %>%
#   mutate(rel_abundance = n / sum(n)) %>%
#   arrange(desc(rel_abundance)) %>%
#   rename(product = products, expected = rel_abundance)
# # observados 
# cases <- cases %>%
#   left_join(rep_bgcs, by = c("Bgcs"="id")) #asignar productos a gccs de la red
# obs <- cases %>%
#   count(rep_bgc) %>%
#   mutate(rel_abundance = n / sum(n)) %>%
#   arrange(desc(rel_abundance)) %>%
#   rename(product = rep_bgc, observed = rel_abundance)
# 
# enrichment <- obs %>%
#   left_join(esp, by = "product") %>%
#   mutate(log2_enrichment = log2(observed / expected)) %>%
#   arrange(desc(log2_enrichment)) 
# 
# # grafica 
# enrichment_plot <- enrichment %>%
#   ggplot(aes(
#     x = log2_enrichment,
#     y = reorder(product, log2_enrichment),
#   )) +
#   geom_col() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   theme_minimal() +
#   labs(
#     x = expression(Log[2](Observed/Expected)),
#     y = "Biosynthetic product"
#   )

### Save tables & graphs ----
write.csv(nodes, paste0(opt$outdir, "nodes_mb.csv"), row.names = FALSE)
write.csv(edges, paste0(opt$outdir, "edges_mb.csv"), row.names = FALSE)
# ggsave(filename = paste0(opt$outdir, "enrichment.png"), plot = enrichment_plot, width = 20, height = 10, units = "cm")

#------------------------------------------
#### MAG->BGC->MAG ####
#------------------------------------------

### Production x Interactions ----
# List of the BGCs and the MAGs that produce them (and their respective groups GCCs and mOTUs)
connection <- meta_bgcs %>%
  select(BGC, Genome, !!bgc_group) %>%
  left_join(meta_mags %>% select(Genome, !!mag_lineage), by = "Genome") %>%
  rename(bgc_id = BGC, gcc_id = !!bgc_group, motu_id = !!mag_lineage) %>%
  distinct()
gcc_in_cases <- unique(cases$Bgcs) # the GCCs that are in the network

# production edges (MAG -> BGC)
production_edges <- connection %>%
  filter(!is.na(motu_id), !is.na(gcc_id)) %>% # keep just the groups
  filter(gcc_id %in% gcc_in_cases) %>%   # keep the ones that are in the network
  select(MAG_prod = motu_id, BGC = gcc_id) %>% # to identify the source
  distinct() 

# interaction edges (BGC -> MAG)
interaction_edges <- edges %>% 
  select(BGC = source, MAG_targ = target, weight)  # to identify the target

# paths of MAGp -> BGC -> MAGi (triplets)
paths <- production_edges %>%
  inner_join(interaction_edges, by = "BGC", relationship = "many-to-many") %>%                             
  filter(MAG_prod != MAG_targ)
# these are all the possible combinations of productions x interactions, we need to filter

### FILTER BY SITES ----
# we need to identify the paths that actually happen in at least 1 site
paths_list <- list()

for (i in seq_len(nrow(paths))) {
  magi <- paths$MAG_prod[i]  # for every pair of mags in a path
  magj <- paths$MAG_targ[i]
  
  comb <- recreate_tableMM(magi, magj, mags_by_sites) # create the table of sites
  # check if they co-occur 
  shared_sites <- sum(comb[[magi]] > 0 & comb[[magj]] > 0) 
  if(shared_sites > 0) {
    paths_list[[i]] <- paths[i, ] %>%
      mutate(oc_sites = shared_sites) # save sites for other filtering maybe?
  }
}

real_paths <- bind_rows(paths_list)
# now we reconstruct the edges table 

### Production x Interactions FILT ----
# reconstruct de production and interaction tables and bind them
production_edges <- real_paths %>%
  select(source = MAG_prod, target = BGC) %>%
  distinct() %>%
  mutate(weight = 1, type = "production")

interaction_edges <- real_paths %>%
  select(source = BGC, target = MAG_targ, weight) %>%
  distinct() %>%
  mutate(type = "interaction")

### NODES AND EDGES ----
edges_mbm <-  bind_rows(production_edges, interaction_edges)

nodes_mbm <- bind_rows(
  edges_mbm %>%
    transmute(id = source, type = if_else(type == "production", "MAG", "BGC")),
  edges_mbm %>%
    transmute(id = target, type = if_else(type == "interaction", "MAG", "BGC")),
) %>%
  distinct()

### GRAPH AND DEGREES ----
g_mbm <- graph_from_data_frame(d = edges_mbm, vertices = nodes_mbm, directed = TRUE)
nodes_mbm$degree_in  <- degree(g_mbm, mode = "in")
nodes_mbm$degree_out <- degree(g_mbm, mode = "out")
nodes_mbm$degree_total <- degree(g_mbm, mode = "all")  

### Add representative groups ----
nodes_mbm <- nodes_mbm %>%
  left_join(rep_mags, by = "id") %>%
  left_join(rep_bgcs, by = "id") %>%
  mutate(color_group = if_else(type == "MAG", rep_mag, rep_bgc)) %>%
  select(-rep_bgc, -rep_mag)

### Save tables ----
write.csv(nodes_mbm, paste0(opt$outdir, "nodes_mbm.csv"), row.names = FALSE)
write.csv(edges_mbm, paste0(opt$outdir, "edges_mbm.csv"), row.names = FALSE)

#-------------------
#### MAG-->MAG ####
#-------------------

### NODES AND EDGES ----
# we use the previous filtered paths to create the edges table
edges_mm <- real_paths %>%
  transmute(source = MAG_prod, target = MAG_targ, bgc = BGC, 
            oc_sites = oc_sites, weight)

nodes_mm <- tibble(id = unique(c(edges_mm$source, edges_mm$target)))

### GRAPH AND DEGREES ----
g_mm <- graph_from_data_frame(d = edges_mm, vertices = nodes_mm, directed = TRUE)
nodes_mm$degree_in  <- degree(g_mm, mode = "in")
nodes_mm$degree_out <- degree(g_mm, mode = "out")
nodes_mm$degree_total <- degree(g_mm, mode = "all")  

### Add representative groups ----
nodes_mm <- nodes_mm %>%
  left_join(rep_mags, by = "id") %>%
  mutate(color_group = rep_mag) %>%
  select(-rep_mag)

### Save tables ----
write.csv(nodes_mm, paste0(opt$outdir, "nodes_mm.csv"), row.names = FALSE)
write.csv(edges_mm, paste0(opt$outdir, "edges_mm.csv"), row.names = FALSE)

