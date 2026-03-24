#######################################
#### NETWORKS MAGs-BGCs ###############
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
  make_option(c("-b", "--bgc_groups"), type="character", default="gcc", help="Name of the grou"),
  make_option(c("-i", "--indir"), type="character", help="Input directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-w", "--workdir"), type="character", help = "Working directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

mag_lineage <- opt$microbial_lineage
bgc_group <- opt$bgc_groups

# Loading metadata and cases
meta_mags <- read.csv(file = paste0(opt$indir, 'metadata.csv'), header = TRUE)
meta_bgcs <- read.csv(file = paste0(opt$indir, 'bgcs_metadata.csv'), header = TRUE)
cases <- read.csv(file= paste0(opt$workdir, 'oc_filt.csv'), header = TRUE)

#----------------
#### MAG<-BGC ####
#----------------

### Nodes and Edges ----
nodes <- tibble(
  id = unique(c(cases$Mags, cases$Bgcs)),
  type = ifelse(id %in% cases$Mags, "MAG", "BGC"))
edges <- cases %>%
  rename(source = Bgcs, target = Mags, weight = fdr_pval_o)

### Graph and degree ----
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
nodes$degree <- degree(g)

### Add representative groups ----
rep_bgcs <- meta_bgcs %>%
  group_by(.data[[bgc_group]], products) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  slice_max(order_by = n, n = 1) %>%
  select(id = all_of(bgc_group), rep_bgc = products)

rep_mags <- meta_mags %>%
  group_by(.data[[mag_lineage]], family) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  slice_max(order_by = n, n = 1) %>%
  select(id = all_of(mag_lineage), rep_mag = family)

nodes <- nodes %>%
  left_join(rep_mags, by = "id") %>%
  left_join(rep_bgcs, by = "id") %>%
  mutate(color_group = if_else(type == "MAG", rep_mag, rep_bgc)) %>%
  select(-rep_bgc, -rep_mag)

### Top nodes ----
top_mag_nodes <- nodes %>%
  filter(type == "MAG") %>%
  arrange(desc(degree)) %>%
  slice_head(n = 12)
top_bgc_nodes <- nodes %>%
  filter(type == "BGC") %>% 
  arrange(desc(degree)) %>%
  slice_head(n = 12)
top_mag_groups <- unique(top_mag_nodes$color_group)
top_bgc_groups <- unique(top_bgc_nodes$color_group)

## Colors for top nodes
mag_colors <- substr(
  as.vector(paletteer::paletteer_d(
    palette = "ggthemes::calc",
    n = length(top_mag_groups))), 1, 7)
bgc_colors <- substr(
  as.vector(paletteer::paletteer_d(
    palette = "ggthemes::gdoc",
    n = length(top_bgc_groups))), 1, 7)

nodes <- nodes %>%
  mutate(color = case_when(
    type == "MAG" & color_group %in% top_mag_groups ~ 
      mag_colors[match(color_group, top_mag_groups)],
    type == "BGC" & color_group %in% top_bgc_groups ~ 
      bgc_colors[match(color_group, top_bgc_groups)],
    TRUE ~ "darkgrey"))


### Crear red en Cytoscape ----
# cytoscapePing()
# createNetworkFromDataFrames(
#  nodes = nodes,
#  edges = edges,
#  title = "oc mOTUs-GCCs",
#  collection = "Interacciones MAGs-BGCs")

### Sub y sobre abundancia productos BGCs  ---

# esperados (dataset global)
esp <- meta_bgcs %>%
  count(products) %>%
  mutate(rel_abundance = n / sum(n)) %>%
  arrange(desc(rel_abundance)) %>%
  rename(product = products, expected = rel_abundance)
# observados 
cases <- cases %>%
  left_join(rep_bgcs, by = c("Bgcs"="id")) #asignar productos a gccs de la red
obs <- cases %>%
  count(rep_bgc) %>%
  mutate(rel_abundance = n / sum(n)) %>%
  arrange(desc(rel_abundance)) %>%
  rename(product = rep_bgc, observed = rel_abundance)

enrichment <- obs %>%
  left_join(esp, by = "product") %>%
  mutate(log2_enrichment = log2(observed / expected)) %>%
  arrange(desc(log2_enrichment)) 

# grafica 
enrichment_plot <- enrichment %>%
  ggplot(aes(
    x = log2_enrichment,
    y = reorder(product, log2_enrichment),
  )) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = expression(Log[2](Observed/Expected)),
    y = "Biosynthetic product"
  )

### Save data & graphs ----
write.csv(nodes, paste0(opt$outdir, "nodes_mb.csv"), row.names = FALSE)
write.csv(edges, paste0(opt$outdir, "edges_mb.csv"), row.names = FALSE)
ggsave(filename = paste0(opt$outdir, "enrichment.png"), plot = enrichment_plot, width = 20, height = 10, units = "cm")

#-----------------------
#### MAG->BGC->MAG ####
#-----------------------

# Connect 
bgc_to_motu <- meta_bgcs %>%
  select(BGC, Genome, !!bgc_group) %>%
  left_join(meta_mags %>% select(Genome, !!mag_lineage), by = "Genome") %>%
  rename(bgc_id = BGC, gcc_id = !!bgc_group, motu_id = !!mag_lineage) %>%
  distinct()
gcc_in_cases <- unique(cases$Bgcs) # keep only the GCCs that are in the network


### Edges ----
# production edges (MAG -> BGC)
production_edges <- bgc_to_motu %>%
  filter(!is.na(motu_id), !is.na(gcc_id)) %>%
  filter(gcc_id %in% gcc_in_cases) %>%
  select(source = motu_id, target = gcc_id) %>%
  distinct() %>%
  mutate(weight = 1, type = "production")

# interaction edges (BGC -> MAG)
interaction_edges <- edges %>%
  mutate(type = "interaction")

# Combine edges
edges_mbm <- bind_rows(production_edges, interaction_edges)


### Nodes ----
nodes_mbm <- tibble(
  id = unique(c(edges_mbm$source, edges_mbm$target))) %>%
  mutate(
    type = case_when(
      id %in% meta_mags[[mag_lineage]] ~ "MAG",
      id %in% meta_bgcs[[bgc_group]] ~ "BGC",
      TRUE ~ "other"))

### Graph and degrees ----
g_mbm <- graph_from_data_frame(d = edges_mbm, vertices = nodes_mbm, directed = TRUE)
nodes_mbm$degree_in  <- degree(g_mbm, mode = "in")
nodes_mbm$degree_out <- degree(g_mbm, mode = "out")
nodes_mbm$degree_total <- degree(g_mbm, mode = "all")

### Representative groups & colors ----
nodes_mbm <- nodes_mbm %>%
  left_join(rep_mags, by = "id") %>%
  left_join(rep_bgcs, by = "id") %>%
  mutate(color_group = if_else(type == "MAG", rep_mag, rep_bgc)) %>%
  select(-rep_bgc, -rep_mag)
nodes_mbm <- nodes_mbm %>%
  mutate(color = case_when(
    type == "MAG" & color_group %in% top_mag_groups ~ 
      mag_colors[match(color_group, top_mag_groups)],
    type == "BGC" & color_group %in% top_bgc_groups ~ 
      bgc_colors[match(color_group, top_bgc_groups)],
    TRUE ~ "darkgrey"))

### Save data ----
write.csv(nodes_mbm, paste0(opt$outdir, "nodes_mbm.csv"), row.names = FALSE)
write.csv(edges_mbm, paste0(opt$outdir, "edges_mbm.csv"), row.names = FALSE)

#----------------
#### MAG-MAG ####
#----------------

### Edges ----
# interaction: MAG <- BGC
inter_mm <- interaction_edges %>%
  select(MAG = target, BGC = source, weight = weight)
# production: MAG -> BGC
prod_mm <- production_edges %>%
  select(MAG = source, BGC = target)
# keep only MAGs
edges_mm <- inter_mm %>%
  inner_join(prod_mm, by = "BGC", relationship = "many-to-many") %>%
  filter(MAG.x != MAG.y) %>%
  transmute(source = MAG.y, target = MAG.x, weight = weight, bgc = BGC) 

### obtain nodes from edges ----
nodes_mm <- tibble(id = unique(c(edges_mm$source, edges_mm$target)))
nodes_mm <- nodes_mm %>%
  left_join(meta_mags %>%
      select(id = all_of(mag_lineage), family) %>%
      distinct(),
    by = "id")
nodes_mm <- nodes_mm %>%
  distinct(id, .keep_all = TRUE)

### add degrees ----
g_mm <- graph_from_data_frame(d = edges_mm, vertices = nodes_mm, directed = TRUE)
nodes_mm$degree_in  <- degree(g_mm, mode = "in")
nodes_mm$degree_out <- degree(g_mm, mode = "out")
nodes_mm$degree_total <- degree(g_mm, mode = "all")

### top nodes and colors ----
top_nodes <- nodes_mm %>% 
  arrange(desc(degree_total)) %>%
  slice_head(n = 12)
top_families <- unique(top_nodes$family)
motu_colors <- substr(
  as.vector(paletteer::paletteer_d(
    palette = "ggthemes::calc",
    n = length(top_families))), 1, 7)

nodes_mm <- nodes_mm %>%
  mutate(color = case_when(
    family %in% top_families ~ 
      motu_colors[match(family, top_families)],
    TRUE ~ "darkgrey"))

### save data ----
write.csv(nodes_mm, paste0(opt$outdir, "nodes_mm.csv"), row.names = FALSE)
write.csv(edges_mm, paste0(opt$outdir, "edges_mm.csv"), row.names = FALSE)



