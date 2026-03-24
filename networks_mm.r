#######################################
#### NETWORKS MAGs-MAGs ###############
## Andrea Zermeño Díaz ################
# march-2026 ##########################

suppressPackageStartupMessages(library(tidyverse))
library(ggplot2)
library(optparse)
library(igraph)
library(paletteer)

# Args
option_list <- list(
  make_option(c("-m", "--microbial_lineage"), type="character", default="mOTUs_Species_Cluster", help="Name of the microbial lienage"),
  make_option(c("-i", "--indir"), type="character", help="Input directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-w", "--workdir"), type="character", help = "Working directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

mag_lineage <- opt$microbial_lineage

# Loading metadata and cases
meta_mags <- read.csv(file = paste0(opt$indir, 'metadata.csv'), header = TRUE)
cases <- read.csv(file= paste0(opt$workdir, 'oc_filt.csv'), header = TRUE)


### Nodes and Edges ----

edges <- cases %>%
  transmute(
    source = MAGi,
    target = MAGj,
    weight = fdr_pval_o
  )
nodes <- tibble(
  id = unique(c(edges$source, edges$target)))

### Graph and degree ----
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
nodes$degree <- degree(g)

### Add representative groups ----
rep_mags <- meta_mags %>%
  group_by(.data[[mag_lineage]], family) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  slice_max(order_by = n, n = 1) %>%
  select(id = all_of(mag_lineage), rep_mag = family)

nodes <- nodes %>%
  left_join(rep_mags, by = "id") %>%
  mutate(color_group = rep_mag) %>%
  select(-rep_mag)

### Top nodes ----
top_mag_nodes <- nodes %>%
  arrange(desc(degree)) %>%
  slice_head(n = 12)
top_mag_groups <- unique(top_mag_nodes$color_group)

## Colors for top nodes
mag_colors <- substr(
  as.vector(paletteer::paletteer_d(
    palette = "ggthemes::calc",
    n = length(top_mag_groups))), 1, 7)

nodes <- nodes %>%
  mutate(color = case_when(
    color_group %in% top_mag_groups ~ 
      mag_colors[match(color_group, top_mag_groups)],
    TRUE ~ "darkgrey"))

### Crear red en Cytoscape ----
# cytoscapePing()
# createNetworkFromDataFrames(
#  nodes = nodes,
#  edges = edges,
#  title = "oc mOTUs-GCCs",
#  collection = "Interacciones MAGs-BGCs")

### Sub y sobre abundancia de familias  ---

# esperados (dataset global)
# esp <- meta_mags %>%
#   count(family) %>%
#   mutate(rel_abundance = n / sum(n)) %>%
#   arrange(desc(rel_abundance)) %>%
#   rename(family = family, expected = rel_abundance)
# 
# cases <- cases %>%
#   left_join(rep_mags, by = c("Mags"="id"))
# # observados
# obs <- cases %>%
#   count(rep_mag) %>%
#   mutate(rel_abundance = n / sum(n)) %>%
#   arrange(desc(rel_abundance)) %>%
#   rename(family = rep_mag, observed = rel_abundance)
# 
# enrichment <- obs %>%
#   left_join(esp, by = "family") %>%
#   mutate(log2_enrichment = log2(observed / expected)) %>%
#   arrange(desc(log2_enrichment)) 
# 
# # grafica 
# enrichment_plot <- enrichment %>%
#   ggplot(aes(
#     x = log2_enrichment,
#     y = reorder(family, log2_enrichment),
#   )) +
#   geom_col() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   theme_minimal() +
#   labs(
#     x = expression(Log[2](Observed/Expected)),
#     y = "Family"
#   )

### Save data & graphs ----
write.csv(nodes, paste0(opt$outdir, "nodes_mm.csv"), row.names = FALSE)
write.csv(edges, paste0(opt$outdir, "edges_mm.csv"), row.names = FALSE)
# ggsave(filename = paste0(opt$outdir, "enrichment.png"), plot = enrichment_plot, width = 20, height = 10, units = "cm")






