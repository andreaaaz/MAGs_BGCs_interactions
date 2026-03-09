################################################
### Identifying interactions MAG-MAG ##########
## Andrea Zermeño Díaz #########################
# march-2026 ################################

# libraries
suppressPackageStartupMessages(library(tidyverse))
library(optparse)

########### Functions #################

prep_mags <- function(meta_mags, mags){
  mags <- sym(mags)
  # data prep
  num_meta <- nrow(meta_mags)
  meta_mags <- meta_mags %>% filter(!is.na(!!mags)) # omit NAs
  fil_meta <- num_meta - nrow(meta_mags)
  message("\n NOTE:", fil_meta, " MAGs were discarded due to missing ", mag_lineage," classification")
  
  # Count the number of MAGs per site and phylogenetic groups
  mags_by_sites <- meta_mags %>%
    count(sites, !!mags) %>%
    pivot_wider(names_from = !!mags, values_from = n, values_fill = 0)
  
  return(mags_by_sites)
}

recreate_table <- function(mag1, mag2, m_by_sites) {
  # bind col1 and col2
  table_comb <- m_by_sites[, c("sites", mag1, mag2), drop=FALSE]

  # filt rows were both are 0
  table_comb <- table_comb[!(table_comb[[mag1]] == 0 & table_comb[[mag2]] == 0), ]
  
  return(table_comb)
}

#### DATA LOAD ####

# Args
option_list <- list(
  make_option(c("-m", "--microbial_lineage"), type="character", default="mOTUs_Species_Cluster", help="Name of the microbial lienage"),
  make_option(c("-s", "--minimum_sites"), type="numeric", default=10, help="Minimum number of sites where a group is present"),
  make_option(c("-i", "--indir"), type="character", help="Working directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-t", "--temp"), type="character", default="global", help="Range of temperature (max, mid and min)")
)
opt <- parse_args(OptionParser(option_list=option_list))

mag_lineage <- opt$microbial_lineage
min_sites <- opt$minimum_sites
temp_r <- opt$temp

# Run script
# Rscript -m mOTUs_Species_Cluster -s 5 -i /mnt/atgc-d3/sur/users/azermeno/exp/MAGs_BGCs_interactions/
# -o /mnt/atgc-d3/sur/users/azermeno/exp/2025-interacions/

message("\n Preparing input, please wait ...")

meta_mags <- read.csv(file = paste0(opt$workdir, 'metadata.csv'), header = TRUE)
meta_sites <- read.csv(file = paste0(opt$workdir, 'meta_sites.csv'), header = TRUE)

##### TEMPERATURE ######

# definir el rango
if (temp_r == "low") temp_r <- c(-2, 9)
if (temp_r == "mid") temp_r <- c(10, 20)
if (temp_r == "high") temp_r <- c(21, 35)
if (temp_r != "global") {
  meta_sites <- meta_sites %>%   # filtrar sitios por temperatura
    filter(between(temperature_..C., temp_r[1], temp_r[2]))
}
# filtrar MAGs por sitios
meta_mags <- meta_mags %>%           
  semi_join(meta_sites, by = "sites")

# Count how many microbial lineages and BGC groups are per site
mags_by_sites <- prep_mags(meta_mags, mag_lineage)

############### Workflow ##############################

# tabla donde se va a guardar informacion de los patrones
cases_list <- list()
total_sites <- length(mags_by_sites)
# para imprimir avance
start_time <- Sys.time()
# para recorrer las columnas de mags_by_sites
mag_cols <- colnames(mags_by_sites)
mag_cols <- mag_cols[mag_cols != "sites"]
n_mags <- length(mag_cols)
# esto generara 'm x m' tablas de 806 renglones (sitios definidos por station_depth)
# donde 'm' es el numero de microbial lienages
for (i in 1:(n_mags - 1)) {
  
  magi <- mag_cols[i]
  
  if (i %% 100 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    message("\n- Microbial lineages processed:", i)
    message("- Time: ", round(elapsed, 2), " mins ...")
  }
  
  for (j in (i + 1):n_mags) {
    
    magj <- mag_cols[j]
    
    comb <- recreate_table(magi, magj, mags_by_sites)
    
    if (nrow(comb) == 0) next 
    
    magi_sites <- sum(comb[[magi]] > 0) 
    magj_sites <- sum(comb[[magj]] > 0)
    
    if (magi_sites < min_sites || magj_sites < min_sites) next
    
    q <- magi_sites / total_sites
    p <- magj_sites / total_sites
    
    pi_e <- p - 2*p*q + q
    pi_o <- p * q
    
    ex_sites <- sum((comb[[magi]] == 0 & comb[[magj]] > 0) |
                      (comb[[magi]] > 0 & comb[[magj]] == 0))
    oc_sites <- sum((comb[[magi]] > 0) & 
                      (comb[[magj]] > 0))
    
    cases_list[[length(cases_list) + 1]] <- list(
      MAGi = magi,
      MAGj = magj,
      mag1_sites = magi_sites,
      mag2_sites = magj_sites,
      ex_sites = ex_sites,
      oc_sites = oc_sites,
      q = q,
      p = p,
      pi_exclusion = pi_e,
      pi_occurrence = pi_o,
      pvalue_e = 1 - pbinom(ex_sites, total_sites, pi_e),
      pvalue_o = 1 - pbinom(oc_sites, total_sites, pi_o)
    )
  }
}
message("\n DONE :)")

# list to data frame
message("\n Preparing output, please wait ...")
cases <- bind_rows(cases_list)

# Correcting multiple testing FDR

cases <- cases %>%
  mutate(
    fdr_pval_e = p.adjust(pvalue_e, method = "BH"),
    fdr_pval_o = p.adjust(pvalue_o, method = "BH")
  )

occurrence <- cases %>%
  filter(fdr_pval_o <= 0.05)
exclusion <- cases %>%
  filter(fdr_pval_e <= 0.05)

# save the produced tables
write.csv(cases, file = paste0(opt$outdir, 'all_cases.csv'), row.names = FALSE)
write.csv(occurrence, file = paste0(opt$outdir, 'oc_filt.csv'), row.names = FALSE)
write.csv(exclusion, file = paste0(opt$outdir, 'ex_filt.csv'), row.names = FALSE)

message("\n Output saved")




