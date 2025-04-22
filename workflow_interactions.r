################################################
### Identifying interactions workflow ##########
## Andrea Zermeño Díaz #########################
# marzo-2025 ###################################

# libraries
library(tidyverse)
library(dplyr)

# data load
workdir <- "/mnt/atgc-d3/sur/users/azermeno/data/2025-interactions/"
outdir <- "/mnt/atgc-d3/sur/users/azermeno/exp/2025-interacions/"

meta_mags <- read.csv(file = paste0(workdir, 'metadata.csv'), header = TRUE)
meta_bgcs <- read.csv(file = paste0(workdir, 'bgcs_metadata.csv'), header = TRUE)

# data prep
# We are grouping by mOTUS so we need to omit the NAs / For other categories skip this 
meta_mags <- meta_mags %>% filter(!is.na(mOTUs.Species.Cluster))
meta_bgcs <- meta_bgcs %>%
  filter(Genome %in% meta_mags$Genome)

# number of mOTUs
#motus <- meta_mags %>% 
#  distinct(mOTUs.Species.Cluster)

# Count the number of MAGs per site and mOTUs
motus_by_sites <- meta_mags %>%
  count(station, mOTUs.Species.Cluster) %>%
  pivot_wider(names_from = mOTUs.Species.Cluster, values_from = n, values_fill = 0)

# bind the column station to each BGC using the column Genome ass a reference to match it 
bgcs_sites <- meta_bgcs %>%
  left_join(meta_mags %>% select(Genome, station), by = "Genome")

# number of GCFs 
#gcfs <- bgcs_sites %>% 
#  distinct(gcf)

# count the number of BGCs per site and GCF
gcfs_by_sites <- bgcs_sites %>%
  count(station, gcf) %>%
  pivot_wider(names_from = gcf, values_from = n, values_fill = 0)


############### Workflow ##############################

# Function
recreate_table <- function(mag, bgc, motus_by_sites, gcfs_by_sites) {
  
  # select the station column and the current column for both tables
  table1 <- motus_by_sites[, c("station", mag), drop = FALSE]
  table2 <- gcfs_by_sites[, c("station", bgc), drop = FALSE]
  
  # bind col1 and col2 by station (full_join porque sino se pierden interacciones)
  table_comb <- full_join(table1, table2, by = "station")
  
  # convertir a 0 los NAs generados por el join
  table_comb[[mag]] <- ifelse(is.na(table_comb[[mag]]), 0, table_comb[[mag]])
  table_comb[[bgc]] <- ifelse(is.na(table_comb[[bgc]]), 0, table_comb[[bgc]])
  
  # Eliminar filas donde ambos valores sean 0
  table_comb <- table_comb[!(table_comb[[mag]] == 0 & table_comb[[bgc]] == 0), ]
  
  return(table_comb)
}

# tablas donde se van a guardar informacion de los patrones
co_occur_list <- list()
co_exclu_list <- list()

counter <- 0
start_time <- Sys.time()

for (col1 in colnames(motus_by_sites)) {
  if (col1 == "station") next
  
  # Imprimir progreso
  counter <- counter + 1
  if (counter %% 10 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    # percentage <- round((counter / ncol(motus_by_site)) * 100, 2)
    
    message("\n- MAGs procesados:", counter)
    message("- Tiempo de proceso: ", round(elapsed, 2), " mins")
    
  }
  
  for (col2 in colnames(gcfs_by_sites)) {
    if (col2 == "station") next
    
    # llamar a la funcon
    temp3 <- recreate_table(col1, col2, motus_by_sites, gcfs_by_sites)
    
    if (nrow(temp3) == 0) next 
    
    # co-ocurrencia
    if (all((temp3[[col1]] > 0) == (temp3[[col2]] > 0))) {
      co_occur_list[[length(co_occur_list) + 1]] <- list(
        mags = col1, 
        bgcs = col2, 
        total_sites = nrow(temp3)
      )
    }
    
    # co-exclusión
    if (all((temp3[[col1]] == 0 & temp3[[col2]] > 0) | (temp3[[col1]] > 0 & temp3[[col2]] == 0))) {
      co_exclu_list[[length(co_exclu_list) + 1]] <- list(
        mags = col1, 
        bgcs = col2, 
        total_sites = nrow(temp3), 
        mag_sites = sum(temp3[[col1]] > 0), 
        bgc_sites = sum(temp3[[col2]] > 0)
      )
    }
  }
}

# Convertir listas a data frames
co_occur <- bind_rows(co_occur_list)
co_exclu <- bind_rows(co_exclu_list)



# Save the produced tables
write.csv(co_exclu, file = paste0(outdir, 'co_exclussion.csv'), row.names = FALSE)
write.csv(co_occur, file = paste0(outdir, 'co_ocurrence.csv'), row.names = FALSE)