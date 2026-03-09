################################################
#### NETWORKS #######################
## Andrea Zermeño Díaz #########################
# march-2026 #################################

# libraries
library(tidyverse)
library(ggplot2)
library(igraph)
library(RCy3)
library(paletteer)


# --- mOTUs - GCC co-ocurrencia

oc <- motu_gcc$occurrence

### 1. Colores por grupo representativo ----
rep_bgcs <- meta_bgcs %>%
  group_by(gcc, products) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  slice_max(order_by = n, n = 1) %>%
  select(gcc, rep_bgc = products) %>%
  rename(id = gcc)
rep_mags <- meta_mags %>%
  group_by(mOTUs_Species_Cluster, family) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  slice_max(order_by = n, n = 1) %>%
  select(mOTUs_Species_Cluster, rep_mag = family) %>%
  rename(id = mOTUs_Species_Cluster)

### 2. Crear nodos y aristas ----
nodes <- data.frame(
  id = unique(c(oc$Mags, oc$Bgcs)),
  type = c(rep("MAG", length(unique(oc$Mags))),
           rep("BGC", length(unique(oc$Bgcs))))
)
edges <- oc %>%
  rename(source = Mags, target = Bgcs, weight = fdr_pval_o)

### 3. Crear grafo y calcular grados ----
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
nodes$degree <- degree(g)

### 4. Añadir metadatos ----
nodes <- nodes %>%
  left_join(rep_mags, by = "id") %>%
  left_join(rep_bgcs, by = "id") %>%
  mutate(color_group = if_else(type == "MAG", rep_mag, rep_bgc))

### 5. Seleccionar top grupos por grado medio ----
top_mag_nodes <- nodes %>%
  filter(type == "MAG") %>%
  arrange(desc(degree)) %>%
  slice_head(n = 13)
top_bgc_nodes <- nodes %>%
  filter(type == "BGC") %>% 
  arrange(desc(degree)) %>%
  slice_head(n = 14)
top_mag_groups <- unique(top_mag_nodes$color_group)
top_bgc_groups <- unique(top_bgc_nodes$color_group)

### 6. Crear red en Cytoscape ----
cytoscapePing()
createNetworkFromDataFrames(
  nodes = nodes,
  edges = edges,
  title = "oc mOTUs-GCCs",
  collection = "Interacciones MAGs-BGCs"
)

### 7. Asignar forma de nodos ----
setNodeShapeMapping("type", c("MAG", "BGC"), c("ELLIPSE", "DIAMOND"))

### 8. Paletas de colores ----
mag_colors <- substr(
  as.vector(paletteer::paletteer_d(
    palette = "ggthemes::calc",
    n = length(top_mag_groups)
  )),
  1, 7
)
bgc_colors <- substr(
  as.vector(paletteer::paletteer_d(
    palette = "ggthemes::gdoc",
    n = length(top_bgc_groups)
  )),
  1, 7
)

### 9. Asignar color a nodos ----
nodes <- nodes %>%
  mutate(color = case_when(
    type == "MAG" & color_group %in% top_mag_groups ~ 
      mag_colors[match(color_group, top_mag_groups)],
    type == "BGC" & color_group %in% top_bgc_groups ~ 
      bgc_colors[match(color_group, top_bgc_groups)],
    TRUE ~ "darkgrey"
  ))

### 10. Aplicar colores en Cytoscape ----
setNodeColorBypass(
  node.names = nodes$id,
  new.colors = nodes$color
)

# guia leyenda
legend_df <- tibble(
  Tipo = c(rep("MAG", length(top_mag_groups)),
           rep("BGC", length(top_bgc_groups))),
  Grupo = c(top_mag_groups, top_bgc_groups),
  Color = c(mag_colors, bgc_colors)
)


# ---- Exclusion mOTUs-GCCs -----

ex <- motu_gcc$exclusion

### 1. Colores por grupo representativo ----
rep_bgcs <- meta_bgcs %>%
  group_by(gcc, products) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  slice_max(order_by = n, n = 1) %>%
  select(gcc, rep_bgc = products) %>%
  rename(id = gcc)
rep_mags <- meta_mags %>%
  group_by(mOTUs_Species_Cluster, family) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  slice_max(order_by = n, n = 1) %>%
  select(mOTUs_Species_Cluster, rep_mag = family) %>%
  rename(id = mOTUs_Species_Cluster)

### 2. Crear nodos y aristas ----
nodes <- data.frame(
  id = unique(c(ex$Mags, ex$Bgcs)),
  type = c(rep("MAG", length(unique(ex$Mags))),
           rep("BGC", length(unique(ex$Bgcs))))
)
edges <- ex %>%
  rename(source = Mags, target = Bgcs, weight = fdr_pval_e)

### 3. Crear grafo y calcular grados ----
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
nodes$degree <- degree(g)

### 4. Añadir metadatos ----
nodes <- nodes %>%
  left_join(rep_mags, by = "id") %>%
  left_join(rep_bgcs, by = "id") %>%
  mutate(color_group = if_else(type == "MAG", rep_mag, rep_bgc))

### 5. Seleccionar top grupos por grado medio ----
top_mag_nodes <- nodes %>%
  filter(type == "MAG") %>%
  arrange(desc(degree)) %>%
  slice_head(n = 2)
top_bgc_nodes <- nodes %>%
  filter(type == "BGC") %>% 
  arrange(desc(degree)) %>%
  slice_head(n = 12)
top_mag_groups <- unique(top_mag_nodes$color_group)
top_bgc_groups <- unique(top_bgc_nodes$color_group)

### 6. Crear red en Cytoscape ----
createNetworkFromDataFrames(
  nodes = nodes,
  edges = edges,
  title = "mOTUs-GCCs",
  collection = "Interacciones MAGs-BGCs"
)

### 7. Asignar forma de nodos ----
setNodeShapeMapping("type", c("MAG", "BGC"), c("ELLIPSE", "DIAMOND"))

### 8. Paletas de colores ----
mag_colors <- substr(
  as.vector(paletteer::paletteer_d(
    palette = "ggthemes::calc",
    n = length(top_mag_groups)
  )),
  1, 7
)
bgc_colors <- substr(
  as.vector(paletteer::paletteer_d(
    palette = "ggthemes::gdoc",
    n = length(top_bgc_groups)
  )),
  1, 7
)

### 9. Asignar color a nodos ----
nodes <- nodes %>%
  mutate(color = case_when(
    type == "MAG" & color_group %in% top_mag_groups ~ 
      mag_colors[match(color_group, top_mag_groups)],
    type == "BGC" & color_group %in% top_bgc_groups ~ 
      bgc_colors[match(color_group, top_bgc_groups)],
    TRUE ~ "darkgrey"
  ))

### 10. Aplicar colores en Cytoscape ----
setNodeColorBypass(
  node.names = nodes$id,
  new.colors = nodes$color
)

# guia leyenda
legend_df <- tibble(
  Tipo = c(rep("MAG", length(top_mag_groups)),
           rep("BGC", length(top_bgc_groups))),
  Grupo = c(top_mag_groups, top_bgc_groups),
  Color = c(mag_colors, bgc_colors)
)

# ---- Sub y sobre abundancia productos BGCs

# esperados (dataset global)
esp <- meta_bgcs %>%
  count(products) %>%
  mutate(rel_abundance = n / sum(n)) %>%
  arrange(desc(rel_abundance)) %>%
  rename(product = products, expected = rel_abundance)

# observados 
oc <- oc %>%
  left_join(rep_bgcs, by = c("Bgcs"="id")) #asignar productos a gccs de la red
obs <- oc %>%
  count(rep_bgc) %>%
  mutate(rel_abundance = n / sum(n)) %>%
  arrange(desc(rel_abundance))%>%
  rename(product = rep_bgc, observed = rel_abundance)
enrichment <- obs %>%
  left_join(esp, by = "product") %>%
  mutate(log2_enrichment = log2(observed / expected)) %>%
  arrange(desc(log2_enrichment)) %>%
  left_join(legend_df, by = c("product" = "Grupo")) %>%
  mutate(
    fill_color = ifelse(is.na(Color), "grey", Color)  
  )

# grafica

enrichment_plot <- enrichment %>%
  ggplot(aes(
    x = log2_enrichment,
    y = reorder(product, log2_enrichment),
    fill = fill_color
  )) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_identity() +
  theme_minimal() +
  labs(
    x = expression(Log[2](Observed/Expected)),
    y = "Biosynthetic product"
  )
ggsave("enrichment.png", plot = enrichment_plot, width = 20, height = 10, units = "cm")




