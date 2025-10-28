################################################
#### ANALYSIS OF RESULTS #######################
## Andrea Zermeño Díaz #########################
# october-2025 #################################

# libraries
library(tidyverse)
library(ggplot2)
library(patchwork)

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
      axis.text.x = element_text(size = 18),  # Tamaño del texto de los números del eje x
      axis.title.x = element_text(size = 18),# Tamaño de la etiqueta del eje y
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
      axis.text.x = element_text(size = 18),  # Tamaño del texto de los números del eje x
      axis.title.x = element_text(size = 18), # Tamaño de la etiqueta del eje y
      panel.background = element_rect(fill = NA, color = NA),
      plot.background = element_rect(fill = NA, color = NA)
    )
  
  list(p_exclusion, p_occurrence)
}


#load data
motu_gcf <- read.csv(file = '~/2025-interacions/motu_gcf/all_cases.csv', header = TRUE)
motu_gcc <- read.csv(file = '~/2025-interacions/motu_gcc/all_cases.csv', header = TRUE)
fam_gcf <- read.csv(file = '~/2025-interacions/fam_gcf/all_cases.csv', header = TRUE)
fam_gcc <- read.csv(file = '~/2025-interacions/fam_gcc/all_cases.csv', header = TRUE)
gen_gcf <- read.csv(file = '~/2025-interacions/gen_gcf/all_cases.csv', header = TRUE)
gen_gcc <- read.csv(file = '~/2025-interacions/gen_gcc/all_cases.csv', header = TRUE)
setwd("~/MAGs_BGCs_interactions")
meta_mags <- read.csv("metadata.csv")
meta_bgcs <- read.csv("bgcs_metadata.csv")

# create the plots
plots <- list(
  plot_pvalues(motu_gcc, 'mOTUs', 'GCCs'),
  plot_pvalues(motu_gcc, 'mOTUs', 'GCCs'),
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

ggsave("/home/azermeno/dist_pvalue.png", plot = p_final, width = 9, height = 4.5, units = "in")


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


#### BUG (la funcion ahora genera una lista lol)
signifs_ex <- tibble(
  mOTUs_gcc = sum(motu_gcc$fdr_pval_e < 0.05)
) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "groups",
    values_to = "num_signif"
  )

signifs_oc <- tibble(
  mOTUs_gcf = sum(motu_gcf$fdr_pval_o < 0.05),
  mOTUs_gcc = sum(motu_gcc$fdr_pval_o < 0.05),
  fam_gcf  = sum(fam_gcf$fdr_pval_o < 0.05),
  fam_gcc  = sum(fam_gcc$fdr_pval_o < 0.05),
  gen_gcf  = sum(gen_gcf$fdr_pval_o < 0.05),
  gen_gcc  = sum(gen_gcc$fdr_pval_o < 0.05)
) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "groups",
    values_to = "num_signif"
  )
bp1 <- ggplot(signifs_ex, aes(x = groups, y = num_signif)) +
  geom_col(fill = "darkblue") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Number of significative co-exclusions (FDR < 0.05)",
    x = "Group type",
    y = "Number of significative cases"
  ) +
  guides(fill = "none")

bp2 <- ggplot(signifs_oc, aes(x = grupos, y = num_signif)) +
  geom_col(fill = "steelblue") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Number of significative co-ocurrences (FDR < 0.05)",
    x = "Group type",
    y = "Number of significative cases"
  ) +
  guides(fill = "none")
bp_final <- (bp1 + bp2) +
  plot_layout(ncol = 2, nrow = 1)
bp_final



### Map of the sites with interaction patterns ###

library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

cases <- read.csv(file = '~/2025-interactions/all_cases.csv', header = TRUE)
meta_mags <- read.csv("~/MAGs_BGCs_interactions/metadata.csv")
example <- recreate_table("meta_mOTU_v25_12843", "gcf_871", mags_by_sites, bgcs_by_sites)

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
    meta_mOTU_v25_12843 > 0 & gcf_871 == 0 ~ "meta_mOTU",
    meta_mOTU_v25_12843 == 0 & gcf_871 > 0 ~ "gcf",
    meta_mOTU_v25_12843 > 0 & gcf_871 > 0  ~ "ambos",
    TRUE                                   ~ "ninguno"
  ))
cols <- c(
  "meta_mOTU" = "#1f77b4",  # azul
  "gcf"       = "#d62728",  # rojo
  "ambos"     = "purple",  # morado
  "ninguno"   = "grey60"    # gris
)

# convertir a objeto espacial para terra
pts <- vect(example, geom = c("longitude", "latitude"), crs = "EPSG:4326")
mundo <- ne_countries(scale = "medium", returnclass = "sf")
oceanos <- ne_download(scale = "medium", type = "ocean", category = "physical", returnclass = "sf")
map <- ggplot() +
  geom_sf(data = oceanos, fill = "#b4d3e2", color = NA) +
  geom_sf(data = mundo, fill = "#5b6684", color = "#5b6684") +
  geom_point(data = example, aes(x = longitude, y = latitude, color = tipo),
             alpha = 0.9, size = 2) +
  scale_color_manual(values = cols, name = "Tipo de interacción") +
  theme_minimal() +
  theme(legend.position = "bottom")
map

### INFO OF THE GROUPS ####

#gcf_871 <- meta_bgcs %>%
#  filter(gcf == "gcf_871")
#mOTU12843 <- meta_mags %>%
#  filter(mOTUs_Species_Cluster == "meta_mOTU_v25_12843")


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
  slice_head(n = 12)
top_bgc_nodes <- nodes %>%
  filter(type == "BGC") %>% 
  arrange(desc(degree)) %>%
  slice_head(n = 13)
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