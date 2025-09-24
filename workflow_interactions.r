################################################
### Identifying interactions workflow ##########
## Andrea Zermeño Díaz #########################
# marzo-2025 ###################################

# libraries
library(tidyverse)
library(optparse)

########### Functions #################

prep_mags <- function(meta_mags, mags){
  mags <- sym(mags)
  # data prep
  meta_mags <- meta_mags %>% filter(!is.na(!!mags)) # omit NAs
  
  # Count the number of MAGs per site and phylogenetic groups
  mags_by_sites <- meta_mags %>%
    count(station, !!mags) %>%
    pivot_wider(names_from = !!mags, values_from = n, values_fill = 0)
  
  return(mags_by_sites)
}

prep_bgcs <- function(meta_bgcs, bgcs){
  bgcs <- sym(bgcs)
  
  # count the number of BGCs per site and GCF
  bgcs_by_sites <- meta_bgcs %>%
    count(station, !!bgcs) %>%
    pivot_wider(names_from = !!bgcs, values_from = n, values_fill = 0)
  
  return(bgcs_by_sites)
}  


recreate_table <- function(mag, bgc, m_by_sites, b_by_sites) {
  
  # select the station column and the current column for both tables
  table1 <- m_by_sites[, c("station", mag), drop = FALSE]
  table2 <- b_by_sites[, c("station", bgc), drop = FALSE]
  
  # bind col1 and col2 by station (full_join porque sino se pierden interacciones)
  table_comb <- full_join(table1, table2, by = "station")
  
  # convertir a 0 los NAs generados por el join
  table_comb[[mag]] <- ifelse(is.na(table_comb[[mag]]), 0, table_comb[[mag]])
  table_comb[[bgc]] <- ifelse(is.na(table_comb[[bgc]]), 0, table_comb[[bgc]])
  
  # filt rows were both are 0
  table_comb <- table_comb[!(table_comb[[mag]] == 0 & table_comb[[bgc]] == 0), ]

  return(table_comb)
}


#### DATA LOAD ####

# Args
option_list <- list(
  make_option(c("-m", "--microbial_lineage"), type="character", default="mOTUs_Species_Cluster", help="Name of the microbial lienage"),
  make_option(c("-b", "--bgc_groups"), type="character", default="gcf", help="Name of the grou"),
  make_option(c("-s", "--minimum_sites"), type="numeric", default=5, help="Minimum number of sites where a group is present"),
  make_option(c("-i", "--workdir"), type="character", help="Working directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

mag_lineage <- opt$microbial_lineage
bgc_group <- opt$bgc_groups
min_sites <- opt$minimum_sites

# Run script
# Rscript -m mOTUs_Species_Cluster -b gcf -s 5 -i /mnt/atgc-d3/sur/users/azermeno/exp/MAGs_BGCs_interactions/
# -o /mnt/atgc-d3/sur/users/azermeno/exp/2025-interacions/

meta_mags <- read.csv(file = paste0(opt$workdir, 'metadata.csv'), header = TRUE)
meta_bgcs <- read.csv(file = paste0(opt$workdir, 'bgcs_metadata.csv'), header = TRUE)


# Count how many microbial lineages and BGC groups are per site
mags_by_sites <- prep_mags(meta_mags, mag_lineage)
bgcs_by_sites<- prep_bgcs(meta_bgcs, bgc_group)


############### Workflow ##############################

# tabla donde se va a guardar informacion de los patrones
cases_list <- list()
total_sites <- length(unique(mags_by_sites$station))
#para imprimir avance
counter <- 0
start_time <- Sys.time()

for (col1 in colnames(mags_by_sites)) {
  if (col1 == "station") next
  
  # Imprimir progreso
  counter <- counter + 1
  if (counter %% 10 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    # percentage <- round((counter / ncol(motus_by_site)) * 100, 2)
    
    message("\n- MAGs procesados:", counter)
    message("- Tiempo de proceso: ", round(elapsed, 2), " mins")
    
  }
  
  for (col2 in colnames(bgcs_by_sites)) { 
    if (col2 == "station") next
    
    # create table of magi and bgcj
    temp3 <- recreate_table(col1, col2, mags_by_sites, bgcs_by_sites)
    
    if (nrow(temp3) == 0) next 
    
    mag_sites <- sum(temp3[[col1]] > 0) 
    bgc_sites <- sum(temp3[[col2]] > 0)
    if (mag_sites < min_sites || bgc_sites < min_sites) next # filtrar MAGs y BGCs que aparezcan en más de 5 sitios 
    
    q <- mag_sites / total_sites
    p <- bgc_sites / total_sites
    pi_e <- p - 2*p*q + q
    pi_o <- p * q
    ex_sites <- sum((temp3[[col1]] == 0 & temp3[[col2]] > 0) |  # co-exclusion
                     (temp3[[col1]] > 0 & temp3[[col2]] == 0))
    oc_sites <- sum((temp3[[col1]] > 0) & (temp3[[col2]] > 0))  # co-occurrence
    
    # save all the combinations of MAGi and BGCj and other data
    cases_list[[length(cases_list) + 1]] <- list(
      # names  
      Mags = col1, 
      Bgcs = col2,
      # sites of mags and bgcs  
      mags_sites = mag_sites, 
      bgcs_sites = bgc_sites,
      # patterns 
      ex_sites = ex_sites,
      oc_sites = oc_sites,
      # probs
      q = q,
      p = p,
      pi_exclusion = pi_e,
      pi_occurrence = pi_o,
      # p-values
      pvalue_e = 1-pbinom(ex_sites, total_sites, pi_e),
      pvalue_o = 1-pbinom(oc_sites, total_sites, pi_o)
      )
  }
}

# list to data frame
cases <- bind_rows(cases_list)

# filter interactions where the BGC is in the Genome
meta_bgcs <- meta_bgcs %>%
  left_join(meta_mags %>% select(Genome, all_of(mag_lineage)), by = "Genome") # add lineage to bgc table by genome

cases <- cases %>% 
  anti_join(meta_bgcs, by = c("Mags" = mag_lineage, "Bgcs" = bgc_group)) 



# Save the produced tables
write.csv(cases, file = paste0(opt$outdir, 'all_cases.csv'), row.names = FALSE)



