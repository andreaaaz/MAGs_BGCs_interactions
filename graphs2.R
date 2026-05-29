meta_mags <- read.csv("~/metadata/metadata.csv")
meta_bgcs <- read.csv("~/metadata/bgcs_metadata.csv")
meta_sites <- read.csv("~/metadata/meta_sites.csv")

cases <- read.csv("~/mOTUs_Species_Cluster_gcc/global/oc_filt.csv")
source("~/MAGs_BGCs_interactions/functions.R")
cases_perm <- read.csv("~/cases_perm.csv")
perm_results <- readRDS("~/perm_results.rds")
################################################
#### ANALYSIS OF RESULTS #######################
## Andrea Zermeño Díaz #########################
# october-2025 #################################

# libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
library(purrr)


#### SITES ANALYSIS #####

# obtener los sitios unicos 
meta_sites <- meta_mags %>%
  distinct(sites, .keep_all = TRUE) %>%
  select(sites, temperature_..C., oxygen_.µmol.kg., date, longitude, latitude, depth, depth_layer, station)

# distribucion de temperatura en los sitios
ggplot(meta_sites, aes(x = temperature_..C.)) + 
  geom_histogram(aes(fill = after_stat(x)), bins = 70) + 
  scale_fill_gradientn(
    colors = c(
      "#08306B",  # azul profundo
      "#41B6C4",  # turquesa
      "#FFFFBF",  # amarillo claro
      "#FC8D59",  # naranja
      "#B30000"   # rojo oscuro
    ),
    name = "Temperature"
  ) + 
  geom_vline(
    xintercept = c(10, 20),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Frequency"
  ) +
  theme_minimal() 

ggplot(meta_sites, aes(x = temperature_..C., y = depth)) +
  geom_point(alpha = 0.6) +
  scale_y_reverse() +
  labs(
    x = "Temperature (°C)",
    y = "Depth (m)"
  ) +
  theme_minimal()


write.csv(meta_sites, "~/MAGs_BGCs_interactions/meta_sites.csv", row.names = FALSE)

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
meta_sites <- read.csv("~/MAGs_BGCs_interactions/meta_sites.csv")

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





