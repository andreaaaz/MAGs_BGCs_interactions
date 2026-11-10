################################################
#### SITES ANALYSIS #######################
## Andrea Zermeño Díaz #########################
# october-2026 #################################

# libraries
library(tidyr)
library(ggplot2)
library(vegan)      
library(cluster)    
library(gridExtra)  

meta_mags <- read.csv("~/metadata/metadata.csv")
meta_bgcs <- read.csv("~/metadata/bgcs_metadata.csv")
meta_sites <- read.csv("~/metadata/meta_sites.csv")
source("~/MAGs_BGCs_interactions/functions.R")

meta_sites <- meta_sites %>%
  filter(n_mags >= 8,
         depth < 2000,
         !is.na(temperature_..C.))%>%
  arrange(desc(temperature_..C.)) %>%
  mutate(temp_group = case_when(
    row_number() <= 150 ~ "high",
    row_number() <= 300 ~ "medium",
    TRUE ~ "low"))
meta_mags <- meta_mags %>%           
  semi_join(meta_sites, by = "sites")
site_info <- meta_sites %>%
  select(sites, temperature_..C., depth_layer, depth, temp_group) %>%
  distinct()

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
  mutate(across(c(n_mags, n_scaffolds, n_families, n_genus, n_motus, n_bgcs, n_gcfs, n_gccs), ~ replace_na(.x, 0))) 

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

ggplot(meta_sites, aes(x = temperature_..C., y = depth, color = temp_group)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = c("high" = "#B30000", "medium" = "#FC8D59", "low" = "#41B6C4"), name = "Temperature") +
  scale_y_reverse() +
  labs(x = "Temperature (°C)", y = "Depth (m)") +
  theme_minimal()

# PRINCIPAL COORDINATE ANALYSIS ---------------

# Data prep --------------------------------------------------------------------
mags_by_sites <- prep_mags(meta_mags, mag_lineage)
mags_by_sites <- mags_by_sites %>% # convertirla en presencia/ausencia
  mutate(across(-sites, ~ as.integer(.x > 0)))

sites_matrix <- mags_by_sites %>% # hacerla matriz
  column_to_rownames("sites") %>%
  as.matrix()
sites_matrix <- sites_matrix[rowSums(sites_matrix) > 0, , drop = FALSE] # quitar sitios sin mOTUs


# Jaccard distances ------------------------------------------------------------
jaccard_dist <- vegdist(sites_matrix, method = "jaccard", binary = TRUE)
k_pcoa <- nrow(sites_matrix) - 1
pcoa <- cmdscale(jaccard_dist, k = k_pcoa, eig = TRUE)

# PCos & eigenvalues explained --------------------------------------------------
# calcular varianza explicada para cada PC 
eig <- pcoa$eig
positive_eig <- eig[eig > 0]
variance_explained <- 100 * positive_eig / sum(positive_eig)
eig_df <- data.frame(axis = seq_along(variance_explained), variance = variance_explained, cumulative = cumsum(variance_explained))

ggplot(eig_df %>% filter(axis <= 50), aes(x = axis, y = variance)) +
  geom_col() +
  geom_line() +
  geom_point() +
  scale_x_continuous( breaks = 1:50, labels = 1:50, limits = c(0.5, 50.5)) +
  labs(x = "PCoA axis",  y = "Percentage of positive eigenvalues (%)") +
  theme_classic()

pcoa_coords <- pcoa$points
colnames(pcoa_coords) <- paste0("PCoA", seq_len(ncol(pcoa_coords)))
n_coords <- ncol(pcoa_coords)
pcoa_df <- as.data.frame(pcoa_coords) %>%
  rownames_to_column("sites") %>%
  left_join(site_info, by = "sites")

ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = temp_group, shape = depth_layer)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = c("high" = "#B30000", "medium" = "#FC8D59", "low" = "#41B6C4"),
                     name = "Temperature") +
  labs(x = paste0("PCoA1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PCoA2 (", round(variance_explained[2], 2), "%)"),
       title = "Microbial community composition across sites") +
  theme_classic()


# CLUSTERING -------------------------------------------------------------------
# comparing different methods
site_ids <- as.character(pcoa_df$sites)

# HIERARQUICAL CLUSTER -----
hc_jaccard <- hclust(jaccard_dist, method = "average")
clusters_jaccard <- cutree(hc_jaccard, k = 3)

sil_h <- silhouette(clusters_jaccard,jaccard_dist) # silhouette
pcoa_df$cluster_h <- factor(clusters_jaccard[pcoa_df$site]) # assign clusters 


# PAM ---------------------
pam_result <- pam(jaccard_dist, k = 3, diss = TRUE)
clusters_pam <- pam_result$clustering

pam_result$silinfo$avg.width # silhouette
pcoa_df$cluster_pam <- factor(clusters_pam[pcoa_df$site]) # assign clusters 

# rename PAM clusters to match hierarchical labels
pcoa_df$cluster_pam <- recode(
  as.character(pcoa_df$cluster_pam),
  "3" = "1",
  "2" = "3",
  "1" = "2"
)


# K-MEANS ----------------
set.seed(123)
km_result <- kmeans(pcoa_coords,centers = 3,nstart = 50)
pcoa_df$cluster_kmeans <- factor(km_result$cluster[site_ids]) # assign clusters


# COMPARISON BETWEEN METHODS ---------------------------------------------------

# Silueta
sil_kmeans <- silhouette(km_result$cluster,jaccard_dist)
sil_comparison <- data.frame(method = c("Hierarchical", "PAM", "K-means"),
  silhouette = c(mean(sil_h[, "sil_width"]),
                 pam_result$silinfo$avg.width,
                 mean(sil_kmeans[, "sil_width"])))
print(sil_comparison)

# Visualizar clustering en PCoAs 
hier_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = cluster_h)) +
  geom_point(size = 1, alpha = 0.8) +
  theme_classic() +
  labs(title = "PCoA - Hierarchical clustering",color = "Cluster")
pam_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = cluster_pam)) +
  geom_point(size = 1, alpha = 0.8) +
  theme_classic() +
  labs(title = "PCoA - PAM clustering", color = "Cluster")
kmean_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = cluster_kmeans)) +
  geom_point(size = 1, alpha = 0.8) +
  theme_classic() +
  labs(title = paste0("PCoA - K-means (", n_coords, " coordinates)"), color = "Cluster")

grid.arrange(hier_plot, pam_plot, kmean_plot, nrow = 1)

# Ver distribucion de temperaturas por clusters (boxplots)
clusters_long <- pcoa_df %>%
  select(sites, temperature_..C., cluster_h, cluster_pam, cluster_kmeans) %>%
  pivot_longer(cols = starts_with("cluster_"), names_to = "method", values_to = "cluster") %>%
  mutate(method = recode(method,cluster_h = "Hierarchical", cluster_pam = "PAM", cluster_kmeans = "K-means"))

ggplot(clusters_long, aes(x = cluster, y = temperature_..C., fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~method) +
  theme_classic() +
  labs(title = "Temperature distribution across microbial clusters", x = "Cluster", y = "Temperature (°C)") +
  guides(fill = "none")

cluster_counts <- clusters_long %>%
  group_by(method, cluster) %>%
  summarise(n = n_distinct(sites), .groups = "drop") %>%
  mutate(cluster_label = paste0("Cluster ", cluster, "\n(n = ", n, ")"))


# TEMPERATURE CORRELATION WITH PCoA1 -------------------------------------------

cor_res <- cor.test(pcoa_df$PCoA1, pcoa_df$temperature_..C., method = "pearson", exact = FALSE)
lm_res  <- lm(PCoA1 ~ temperature_..C., data = pcoa_df)
print(cor_res)
summary(lm_res)

p_val <- cor_res$p.value
lab <- paste0("r = ", round(cor_res$estimate, 3),
              "\nR² = ", round(summary(lm_res)$r.squared, 3),
              "\np ", ifelse(p_val < 0.001, "< 0.001", paste0("= ", round(p_val, 3))))
library(patchwork)
temp_cols <- c("high" = "#B30000", "medium" = "#FC8D59", "low" = "#41B6C4")
x_lab     <- paste0("PCoA1 (", round(variance_explained[1], 2), "%)")
x_range   <- range(pcoa_df$PCoA1)   # mismo rango en x para ambos

p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = temp_group, shape = depth_layer)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = temp_cols, name = "Temperature") +
  scale_x_continuous(limits = x_range) +
  labs(x = x_lab,
       y = paste0("PCoA2 (", round(variance_explained[2], 2), "%)"),
       title = "Microbial community composition across sites") +
  theme_classic()

p_cor <- ggplot(pcoa_df, aes(x = PCoA1, y = temperature_..C.)) +
  geom_point(aes(color = temp_group), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = temp_cols, name = "Temperature") +
  scale_x_continuous(limits = x_range) +
  annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.1, vjust = 1.2, size = 4) +
  labs(x = x_lab, y = "Temperature (°C)") +
  guides(color = "none") +
  theme_classic()
# Juntar: uno sobre otro
(p_pcoa / p_cor) +
  plot_layout(heights = c(2, 1), guides = "collect")






filt_pcoa <- pcoa_df %>%
  filter(PCoA1 > 0)
cor_res <- cor.test(filt_pcoa$PCoA3, filt_pcoa$temperature_..C., method = "pearson", exact = FALSE)
lm_res  <- lm(PCoA3 ~ temperature_..C., data = filt_pcoa)
print(cor_res)
summary(lm_res)

p_val <- cor_res$p.value
lab <- paste0("r = ", round(cor_res$estimate, 3),
              "\nR² = ", round(summary(lm_res)$r.squared, 3),
              "\np ", ifelse(p_val < 0.001, "< 0.001", paste0("= ", round(p_val, 3))))
library(patchwork)
temp_cols <- c("high" = "#B30000", "medium" = "#FC8D59", "low" = "#41B6C4")
x_lab     <- paste0("PCoA1 (", round(variance_explained[1], 2), "%)")
x_range   <- range(filt_pcoa$PCoA1)   # mismo rango en x para ambos

p_pcoa <- ggplot(filt_pcoa, aes(x = PCoA1, y = PCoA2, color = temp_group, shape = depth_layer)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = temp_cols, name = "Temperature") +
  scale_x_continuous(limits = x_range) +
  labs(x = x_lab,
       y = paste0("PCoA2 (", round(variance_explained[2], 2), "%)"),
       title = "Microbial community composition across sites") +
  theme_classic()

p_cor <- ggplot(filt_pcoa, aes(x = PCoA1, y = temperature_..C.)) +
  geom_point(aes(color = temp_group), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = temp_cols, name = "Temperature") +
  scale_x_continuous(limits = x_range) +
  annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.1, vjust = 1.2, size = 4) +
  labs(x = x_lab, y = "Temperature (°C)") +
  guides(color = "none") +
  theme_classic()
# Juntar: uno sobre otro
(p_pcoa / p_cor) +
  plot_layout(heights = c(2, 1), guides = "collect")



ggplot(filt_pcoa, aes(x = PCoA3, y = temperature_..C.)) +
  geom_point(aes(color = temp_group), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = temp_cols, name = "Temperature") +
  scale_x_continuous(limits = x_range) +
  annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.1, vjust = 1.2, size = 4) +
  labs(x = x_lab, y = "Temperature (°C)") +
  guides(color = "none") +
  theme_classic()









pcoa_df <- pcoa_df %>%
  left_join(meta_sites %>% select(sites, n_mags) %>% distinct(), by = "sites")
cor.test(pcoa_df$PCoA1, pcoa_df$n_mags, method = "spearman")

# puedes hacer poa1 vs pcoa 2 pero que el color sea por n_mags?





# Gráfica
ggplot(pcoa_df, aes(x = temperature_..C., y = n_mags, color = temp_group)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = temp_cols, name = "Temperature") +
  theme_classic()

