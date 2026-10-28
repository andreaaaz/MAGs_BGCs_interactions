################################################
#### SITES ANALYSIS #######################
## Andrea Zermeño Díaz #########################
# october-2026 #################################

# libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
library(purrr)

meta_mags <- read.csv("~/metadata/metadata.csv")
meta_bgcs <- read.csv("~/metadata/bgcs_metadata.csv")
meta_sites <- read.csv("~/metadata/meta_sites.csv")
source("~/MAGs_BGCs_interactions/functions.R")

# SEQUENCING QUALITY CONTROL ------------------
# MAGs
qc_mags <- meta_mags %>%
  group_by(sites) %>%
  summarise(
    n_mags = n_distinct(Genome),
    n_scaffolds = sum(X._scaffolds, na.rm = TRUE),
    n_families = n_distinct(family, na.rm = TRUE),
    n_genus = n_distinct(genus, na.rm = TRUE),
    n_motus = n_distinct(mOTUs_Species_Cluster, na.rm = TRUE),
    .groups = "drop"
  )

# BGCs
qc_bgcs <- meta_bgcs %>%
  group_by(sites) %>%
  summarise(
    n_bgcs = n_distinct(BGC),
    n_gcfs = n_distinct(gcf, na.rm = TRUE),
    n_gccs = n_distinct(gcc, na.rm = TRUE),
    .groups = "drop"
  )

# Combine
qc_sites <- qc_mags %>%
  full_join(qc_bgcs, by = "sites") %>%
  mutate(
    across(
      c(n_mags, n_scaffolds, n_families, n_genus, n_motus,
        n_bgcs, n_gcfs, n_gccs),
      ~ replace_na(.x, 0)
    )
  ) 


qc_sites %>% 
  ggplot( aes(x=n_mags)) +
  stat_bin(breaks=seq(0,100,1), fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Number of MAGs in sites") +
  theme_minimal() 

qc_sites <- qc_sites %>%
  full_join(meta_sites, by = "sites")

ggplot(qc_sites, aes(x = temperature_..C., y = n_mags)) +
  geom_point(alpha = 0.7, size = 1) +
  theme_minimal()
ggplot(qc_sites, aes(x = depth, y = n_mags)) +
  geom_point(alpha = 0.7, size = 1) +
  theme_minimal()

# TEMPERATURE ---------------------------------
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
  labs(
    x = "Temperature (°C)",
    y = "Frequency"
  ) +
  theme_minimal() 

ggplot(meta_sites, aes(x = temperature_..C., y = depth, color = temperature_..C.)) +
  geom_point(alpha = 1) +
  scale_color_gradientn(
    colors = c(
      "#08306B",  # azul profundo
      "#41B6C4",  # turquesa
      "#FFFFBF",  # amarillo claro
      "#FC8D59",  # naranja
      "#B30000"   # rojo oscuro
    ),
    name = "Temperature"
  ) +
  scale_y_reverse() +
  labs(
    x = "Temperature (°C)",
    y = "Depth (m)"
  ) +
  theme_minimal()


# PRINCIPAL COORDINATE ANALYSIS ---------------
# meta_sites <- meta_sites %>%
#   filter(n_mags >= qc_sites)
mags_by_sites <- prep_mags(meta_mags, "mOTUs_Species_Cluster")
# convertirla en presencia/ausencia
mags_by_sites <- mags_by_sites %>%
  mutate(across(-sites, ~ as.integer(.x > 0)))


