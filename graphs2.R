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
    geom_histogram(binwidth = 0.05, fill = "darkblue") +
    labs(
      title = paste("Co-exclusion of", mag_group, "and", bgc_group),
      x = "p-value",
      y = "Frequency"
    ) +
    theme_minimal()
  
  p_occurrence <- ggplot(cases, aes(x = pvalue_o)) +
    geom_histogram(binwidth = 0.05, fill = "steelblue") +
    labs(
      title = paste("Co-occurrence of", mag_group, "and", bgc_group),
      x = "p-value",
      y = "Frequency"
    ) +
    theme_minimal()
  
  list(p_exclusion, p_occurrence)
}

#load data
motu_gcf <- read.csv(file = '~/2025-interacions/motu_gcf/all_cases.csv', header = TRUE)
motu_gcc <- read.csv(file = '~/2025-interacions/motu_gcc/all_cases.csv', header = TRUE)
fam_gcf <- read.csv(file = '~/2025-interacions/fam_gcf/all_cases.csv', header = TRUE)
fam_gcc <- read.csv(file = '~/2025-interacions/fam_gcc/all_cases.csv', header = TRUE)
gen_gcf <- read.csv(file = '~/2025-interacions/gen_gcf/all_cases.csv', header = TRUE)
gen_gcc <- read.csv(file = '~/2025-interacions/gen_gcc/all_cases.csv', header = TRUE)

# create the plots
plots <- list(
  plot_pvalues(motu_gcf, 'mOTUs', 'GCFs'),
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

p_final



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
  mOTUs_gcf = sum(motu_gcf$fdr_pval_e < 0.05),
  mOTUs_gcc = sum(motu_gcc$fdr_pval_e < 0.05),
  fam_gcf  = sum(fam_gcf$fdr_pval_e < 0.05),
  fam_gcc  = sum(fam_gcc$fdr_pval_e < 0.05),
  gen_gcf  = sum(gen_gcf$fdr_pval_e < 0.05),
  gen_gcc  = sum(gen_gcc$fdr_pval_e < 0.05)
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

# Red no dirigida
# Cada nodo es un grupo de MAGs/BGCs
# El color de la linea es si la interaccion es de exclusion o ocurrencia

library(ggraph)
library(igraph)
library(RCy3)
cytoscapePing()

cytograph <- function(df, name, rep_mags, rep_bgcs) {
  nodes <- data.frame(
    id = unique(c(df$Mags, df$Bgcs)),
    type = c(rep("MAG", length(unique(df$Mags))),
             rep("BGC", length(unique(df$Bgcs))))
  )
  edges <- df %>%
    rename(source = Mags, target = Bgcs, weight = fdr_pval_o)
  g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  nodes$degree <- degree(g)
  
  nodes <- left_join(nodes, rep_mags, by = mOTUs_Species_Cluster)
  nodes <- left_join(nodes, rep_bgcs, by = gcc)
  
  createNetworkFromDataFrames(nodes = nodes, edges = edges,
                              title = name, collection = "Interacions MAGs-BGCs")
  setNodeShapeMapping("type", c("MAG","BGC"), c("ELLIPSE","DIAMOND"))
  
  top_nodes <- nodes %>%
    arrange(desc(degree)) %>%
    slice_head(n = 8)
  
  setNodeColorDefault("#D3D3D3")  
  setNodeColorBypass(node.names = top_nodes$id,
                           new.colors = RColorBrewer::brewer.pal(min(10, 9), "Set1"))
}
  
oc <- motu_gcc$occurrence
network <- cytograph(oc, "mOTUs-GCCs", rep_mags, rep_bgcs)


### colors by group
rep_bgcs <- meta_bgcs %>%
  group_by(gcc, products) %>%          
  summarise(n = n(), .groups = "drop_last") %>%  
  slice_max(order_by = n, n = 1) %>%    
  select(gcc, rep_bgc = products)

rep_mags <- meta_mags %>%
  group_by(mOTUs_Species_Cluster, family) %>%        
  summarise(n = n(), .groups = "drop_last") %>%  
  slice_max(order_by = n, n = 1) %>%    
  select(mOTUs_Species_Cluster, rep_mag = family)










# MAGs = círculos, BGCs = triángulos
setNodeShapeMapping("type", c("MAG","BGC"), c("ELLIPSE","DIAMOND"))
# Tamaño de nodos
setNodeSizeMapping("degree", min(nodes$degree), max(nodes$degree), mapping.type = "c")
# Grosor de aristas
edges$edgeWidth <- -log10(edges$weight)  # o cualquier transformación
setEdgeLineWidthMapping("edgeWidth", min(edges$edgeWidth), max(edges$edgeWidth), 1, 10, mapping.type = "c")

# Color opcional por tipo
setNodeColorMapping("type", c("MAG","BGC"), c("#1f78b4","#33a02c"))
# 
