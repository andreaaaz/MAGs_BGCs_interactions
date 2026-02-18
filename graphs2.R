################################################
#### ANALYSIS OF RESULTS #######################
## Andrea Zermeño Díaz #########################
# october-2025 #################################

# libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
library(purrr)


######### P-VALUES ANALYSIS ########## 
### Statistical signals ###
## p-values distribution 

plot_pvalues <- function(cases, mag_group, bgc_group) {
  p_exclusion <- ggplot(cases, aes(x = pvalue_e)) +
    geom_histogram(binwidth = 0.05, fill = "#044a88") +
    labs(
      title = paste("Co-exclusion of", mag_group, "and", bgc_group),
      x = "p-value",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),  # Tamaño del texto de los números del eje x
      axis.title.x = element_text(size = 12),# Tamaño de la etiqueta del eje y
    )
  
  p_occurrence <- ggplot(cases, aes(x = pvalue_o)) +
    geom_histogram(binwidth = 0.05, fill = "#b4d3e2") +
    labs(
      title = paste("Co-occurrence of", mag_group, "and", bgc_group),
      x = "p-value",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),  # Tamaño del texto de los números del eje x
      axis.title.x = element_text(size = 12), # Tamaño de la etiqueta del eje y
      panel.background = element_rect(fill = NA, color = NA),
      plot.background = element_rect(fill = NA, color = NA)
    )
  
  list(p_exclusion, p_occurrence)
}


#load data
motu_gcf <- read.csv(file = '~/2025-interacions/motu_gcf/all_cases.csv', header = TRUE)
motu_gcc <- read.csv(file = '~/2026-interactions/motu_gcc/all_cases.csv', header = TRUE)
fam_gcf <- read.csv(file = '~/2026-interactions/fam_gcf/all_cases.csv', header = TRUE)
fam_gcc <- read.csv(file = '~/2026-interactions/fam_gcc/all_cases.csv', header = TRUE)
gen_gcf <- read.csv(file = '~/2026-interactions/gen_gcf/all_cases.csv', header = TRUE)
gen_gcc <- read.csv(file = '~/2026-interactions/gen_gcc/all_cases.csv', header = TRUE)
meta_mags <- read.csv("~/MAGs_BGCs_interactions/metadata.csv")
meta_bgcs <- read.csv("~/MAGs_BGCs_interactions/bgcs_metadata.csv")

# create the plots
plots <- list(
  plot_pvalues(motu_gcc, 'mOTUs', 'GCCs'),
  plot_pvalues(motu_gcf, 'mOTUs', 'GCFs'),
  plot_pvalues(fam_gcf,  'families', 'GCFs'),
  plot_pvalues(fam_gcc,  'families', 'GCCs'),
  plot_pvalues(gen_gcf,  'genus', 'GCFs'),
  plot_pvalues(gen_gcc,  'genus', 'GCCs')
)
# flat the list of lists
plots_flat <- do.call(c, plots)

# paste in only one graph
p_final <- wrap_plots(plots_flat, ncol = 4, nrow = 3) +
  plot_annotation(title = "Distribución de p-values")

ggsave("~/2026-interactions/dist.png", plot = p_final, width = 9, height = 4.5, units = "in")


### Multiple testing correction by FDR ### 

# correcting by FDR and filtering the cases with p-value > 0.05
correct <- function(df) {
  df <- df %>%
    mutate(
      fdr_pval_e = p.adjust(pvalue_e, method = "BH"),
      fdr_pval_o = p.adjust(pvalue_o, method = "BH")
    )
  
  list(
    exclusion = df %>% filter(fdr_pval_e < 0.05),
    occurrence = df %>% filter(fdr_pval_o < 0.05)
  )
}

motu_gcf <- correct(motu_gcf)
motu_gcc <- correct(motu_gcc)
fam_gcf <- correct(fam_gcf)
fam_gcc <- correct(fam_gcc)
gen_gcf <- correct(gen_gcf)
gen_gcc <- correct(gen_gcc)


#### Plot laa cantidad de casos para cada combinacion
listas <- list(
  motu_gcf = motu_gcf,
  motu_gcc = motu_gcc,
  fam_gcf  = fam_gcf,
  fam_gcc  = fam_gcc,
  gen_gcf  = gen_gcf,
  gen_gcc  = gen_gcc
)

sum_cases <- imap_dfr(listas, ~ tibble(
  dataset = .y,
  exclusion_n = nrow(.x$exclusion),
  occurrence_n = nrow(.x$occurrence)
))


bp1 <- ggplot(sum_cases, aes(x = dataset, y = exclusion_n)) +
  geom_col(fill = "darkblue") +
  labs(title = "Co-exclusions (FDR < 0.05)", 
       x = "Group type", 
       y = "Number of significative cases") +
  scale_y_continuous(breaks = seq(0, max(sum_cases$exclusion_n), by = 1)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")



bp2 <- ggplot(sum_cases, aes(x = dataset, y = occurrence_n)) +
  geom_col(fill = "steelblue") +
  labs(title = "Co-occurrences (FDR < 0.05)", 
       x = "Group type", 
       y = "Number of significative cases") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

bp_final <- (bp1 + bp2) +
  plot_layout(ncol = 2, nrow = 1)
bp_final



### Map of the sites with interaction patterns ###

library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

oc<- motu_gcc$occurrence
meta_mags <- read.csv("~/MAGs_BGCs_interactions/metadata.csv")
example <- recreate_table("ref_mOTU_v25_05467", "gcc_22", mags_by_sites, bgcs_by_sites)

#sitios con cordenadas
sites <- data.frame(
  station = meta_mags$station,
  latitude = meta_mags$latitude,
  longitude = meta_mags$longitude
)
sites <- sites %>%
  group_by(station) %>%
  summarise(
    longitude = first(longitude),
    latitude = first(latitude)
  )
#unir sitios a la tabla
example <- example %>%
  left_join(sites %>% select(latitude, longitude, station), by = "station")

example <- example %>%
  mutate(tipo = case_when(
    ref_mOTU_v25_05467 > 0 & gcc_22 == 0 ~ "mOTU",
    ref_mOTU_v25_05467 == 0 & gcc_22 > 0 ~ "GCC",
    ref_mOTU_v25_05467 > 0 & gcc_22 > 0  ~ "Co-occurrence",
    TRUE                                   ~ "NA"
  ))
cols <- c(
  "mOTU" = "#009d14",  
  "GCC" = "#990099",  
  "Co-occurrence" = "darkblue", 
  "ninguno"   = "grey60"   
)

# convertir a objeto espacial para terra
pts <- vect(example, geom = c("longitude", "latitude"), crs = "EPSG:4326")
mundo <- ne_countries(scale = "medium", returnclass = "sf")
oceanos <- ne_download(scale = "medium", type = "ocean", category = "physical", returnclass = "sf")
map <- ggplot() +
  geom_sf(data = oceanos, fill = "#b4d3e2", color = NA) +
  geom_sf(data = mundo, fill = "#5b6684", color = "#5b6684") +
  geom_point(data = example, aes(x = longitude, y = latitude, color = tipo),
             alpha = 0.9, size = 1) +
  scale_color_manual(values = cols, name = "Type of interaction") +
  theme_minimal(base_size = 18) +
  theme( plot.background = element_rect(fill = NA, color = NA),   # fondo transparente
         panel.background = element_rect(fill = NA, color = NA),  # panel transparente
         legend.background = element_rect(fill = NA, color = NA), # fondo de leyenda transparente
         legend.key = element_rect(fill = NA, color = NA),        # clave de leyenda transparente
         axis.text = element_text(size = 16),
         axis.title = element_text(size = 18, face = "bold"),
         legend.title = element_text(size = 18, face = "bold"),
         legend.text = element_text(size = 16))
map
ggsave("map.png", plot = map, width = 42, height = 22, units = "cm", bg = "transparent")

### INFO OF THE GROUPS ####

gcc_22 <- meta_bgcs %>%
  filter(gcc == "gcc_22")
mOTU5467 <- meta_mags %>%
  filter(mOTUs_Species_Cluster == "ref_mOTU_v25_05467")


### NETWORKS ###

library(igraph)
library(RCy3)
library(paletteer)
cytoscapePing()

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


### --- PRODUCTION NETWORK
gcc_22 <- gcc_22 %>%
  separate(Taxonomy,
           into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";",
           fill = "right",
           remove = FALSE) %>%
  mutate(across(
    c(domain, phylum, class, order, family, genus, species),
    ~ str_remove(., "^[a-z]__")  # elimina prefijos  d__ y asi
  )) %>%
  mutate(across(
    c(domain, phylum, class, order, family, genus, species),
    ~ na_if(., "")  # conviertir a NA
  )) 
gcc_22 <- gcc_22 %>%
  left_join(meta_mags %>% select(Genome, mOTUs_Species_Cluster), by = "Genome")

mOTU5467 <- mOTU5467 %>%
  mutate(gcc = "gcc_22")


gcc_node <- gcc_22 %>%
  distinct(gcc) %>%
  transmute(
    id = gcc,
    type = "GCC",
    family = NA
  )
producer_nodes <- gcc_22 %>%
  distinct(mOTUs_Species_Cluster, family) %>%
  transmute(
    id = mOTUs_Species_Cluster,
    type = "Producer",
    family = family
  )
eco_node <- mOTU5467 %>%
  transmute(
    id = mOTUs_Species_Cluster,
    type = "interactor",
    family = family
  )
nodes <- bind_rows(gcc_node, producer_nodes, eco_node) %>%
  distinct(id, .keep_all = TRUE)
production_edges <- gcc_22 %>%
  distinct(mOTUs_Species_Cluster, gcc) %>%
  transmute(
    source = mOTUs_Species_Cluster,
    target = gcc,
    interaction = "production"
  )
eco_edges <- mOTU5467 %>%
  transmute(
    source = mOTUs_Species_Cluster,
    target = gcc,
    interaction = "ecological_interaction"
  )

edges <- bind_rows(production_edges, eco_edges)

createNetworkFromDataFrames(
  nodes = nodes,
  edges = edges,
  title = "GCC production and ecological interaction",
  collection = "GCC ecological networks"
)
setNodeShapeMapping(
  "type",
  c("GCC", "Producer", "Ecological_interactor"),
  c("HEXAGON", "ELLIPSE", "TRIANGLE")
)
families <- nodes %>%
  filter(!is.na(family)) %>%
  pull(family) %>%
  unique()

family_colors <- substr(
  as.vector(
    paletteer::paletteer_d(
      palette = "ggthemes::calc",
      n = length(families)
    )
  ),
  1, 7
)
setNodeColorMapping(
  table.column = "family",
  table.column.values = families,
  colors = family_colors
)
setEdgeLineStyleMapping(
  "interaction",
  c("production", "ecological_interaction"),
  c("SOLID", "DASHED")
)

