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

meta_sites <- meta_sites %>%
  filter(n_mags >= 8)
meta_mags <- meta_mags %>%           
  semi_join(meta_sites, by = "sites")

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
    colors = c("#08306B", "#41B6C4", "#FFFFBF", "#FC8D59", "#B30000"),
    name = "Temperature") + 
  labs(x = "Temperature (°C)", y = "Frequency") +
  theme_minimal() 

ggplot(meta_sites, aes(x = temperature_..C., y = depth, color = temperature_..C.)) +
  geom_point(alpha = 1) +
  scale_color_gradientn(
    colors = c("#08306B", "#41B6C4", "#FFFFBF", "#FC8D59", "#B30000"),
    name = "Temperature") +
  scale_y_reverse() +
  labs(x = "Temperature (°C)", y = "Depth (m)") +
  theme_minimal()


# PRINCIPAL COORDINATE ANALYSIS ---------------
# meta_sites <- meta_sites %>%
#   filter(n_mags >= qc_sites)
mags_by_sites <- prep_mags(meta_mags, mag_lineage)
# convertirla en presencia/ausencia
mags_by_sites <- mags_by_sites %>%
  mutate(across(-sites, ~ as.integer(.x > 0)))
# hacerla matrix
sites_matrix <- mags_by_sites %>%
  column_to_rownames("sites") %>%
  as.matrix()
# quitar sitios sin mOTUs
sites_matrix <- sites_matrix[rowSums(sites_matrix) > 0, , drop = FALSE]

# calcular jaccard distance 
jaccard_dist <- vegdist(sites_matrix, method = "jaccard", binary = TRUE)
pcoa <- cmdscale(jaccard_dist, k = 2, eig = TRUE)
pcoa_scores <- as.data.frame(pcoa$points) %>%
  setNames(c("PCoA1", "PCoA2")) %>%
  rownames_to_column("sites") %>%
  left_join(meta_sites %>% select(sites, temperature_..C., depth_layer, depth) %>% distinct(), by = "sites")

# PCs
eig <- pcoa$eig
# proportion of positive eigenvalues explained
positive_eig <- eig[eig > 0]
variance_explained <- 100 * positive_eig / sum(positive_eig)
# Percentage explained by the first two axes
eig_df <- data.frame(axis = seq_along(variance_explained), variance = variance_explained, cumulative = cumsum(variance_explained))
ggplot(eig_df %>% filter(axis <= 50), aes(x = axis, y = variance)) +
  geom_col() +
  geom_line() +
  geom_point() +
  scale_x_continuous( breaks = 1:50, labels = 1:50, limits = c(0.5, 50.5)) +
  labs(x = "PCoA axis",  y = "Percentage of positive eigenvalues (%)") +
  theme_classic()


# plot PCoA
ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2, color = temperature_..C., shape = depth_layer)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_gradientn(
    colors = c("#08306B", "#41B6C4", "#FFFFBF", "#FC8D59", "#B30000"),
    name = "Temperature") +
  labs(x = paste0("PCoA1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PCoA2 (", round(variance_explained[2], 2), "%)"),
       title = "Microbial community composition across sites") +
  theme_classic()


# herarcical clustering
hc_jaccard <- hclust(jaccard_dist, method = "average")
clusters_jaccard <- cutree(hc_jaccard, k = 3)
sil <- silhouette(clusters_jaccard, jaccard_dist)
mean(sil[, "sil_width"])
aggregate(sil[, "sil_width"], by = list(cluster = sil[, "cluster"]), FUN = mean)

pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$site <- rownames(pcoa_df)
pcoa_df$cluster_h <- factor(clusters_jaccard[pcoa_df$site])

pcoa_df <- pcoa_df %>%
  left_join(meta_sites %>% select(sites, temperature_..C., depth_layer), by = c("site" = "sites"))

ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = cluster_h)) +
  geom_point(size = 1, alpha = 0.8) +
  theme_classic() +
  labs(title = "PCoA - Jaccard", color = "Cluster", x = "PCoA 1", y = "PCoA 2")


# PAM
pam_result <- pam(jaccard_dist, k = 3, diss = TRUE)
clusters_pam <- pam_result$clustering
pam_result$silinfo$avg.width
pcoa_df$cluster_pam <- factor(clusters_pam[match(pcoa_df$site, names(clusters_pam))])

ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = cluster_pam)) +
  geom_point(size = 1, alpha = 0.8) +
  theme_classic() +
  labs(title = "PCoA - PAM clustering",
       color = "Cluster")
table(clusters_pam)


pcoa_df$cluster_pam <- recode(
  as.character(pcoa_df$cluster_pam),
  "1" = "3",
  "2" = "1",
  "3" = "2"
)

pcoa_df$cluster_pam <- factor(
  pcoa_df$cluster_pam,
  levels = c("1", "2", "3")
)
