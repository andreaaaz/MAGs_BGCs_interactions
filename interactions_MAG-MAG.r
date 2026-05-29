################################################
### Identifying interactions MAG-MAG ##########
## Andrea Zermeño Díaz #########################
# march-2026 ################################

# libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
library(optparse)

# Args
option_list <- list(
  make_option(c("-m", "--microbial_lineage"), type="character", default="mOTUs_Species_Cluster", help="Name of the microbial lienage"),
  make_option(c("-s", "--minimum_sites"), type="numeric", default=10, help="Minimum number of sites where a group is present"),
  make_option(c("-i", "--indir"), type="character", help="Input directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-w", "--workdir"), type="character", help="Working directory"),
  make_option(c("-t", "--temp"), type="character", default="high", help="Range of temperature (max, mid and min)"),
  make_option(c("-e", "--method"), type="character", default="binomial", help="Method to calculate significance")
)
opt <- parse_args(OptionParser(option_list=option_list))

mag_lineage <- opt$microbial_lineage
min_sites <- opt$minimum_sites
temp_r <- opt$temp
method <- opt$method

# Run script
# Rscript -m mOTUs_Species_Cluster -s 5 -i /mnt/atgc-d3/sur/users/azermeno/exp/MAGs_BGCs_interactions/
# -o /mnt/atgc-d3/sur/users/azermeno/exp/2025-interacions/


#### DATA LOAD ####
message("\n Preparing input, please wait ...")

meta_mags <- read.csv(file = paste0(opt$indir, 'metadata.csv'), header = TRUE)
meta_sites <- read.csv(file = paste0(opt$indir, 'meta_sites.csv'), header = TRUE)
# functions
source(paste0(opt$workdir, "functions.R"))


##### TEMPERATURE ######
temp_range <- NULL

# definir el rango
if (temp_r == "low") temp_range <- c(-2, 9)
if (temp_r == "mid") temp_range <- c(10, 20)
if (temp_r == "high") temp_range <- c(21, 35)
if (temp_r != "global") {
  meta_sites <- meta_sites %>%   # filtrar sitios por temperatura
    filter(between(temperature_..C., temp_range[1], temp_range[2]))
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
    
    res <- binomial_MM(magi, magj, mags_by_sites, min_sites, total_sites)
    
    # choosing method
    if (method == "binomial") {
      res <- binomial_MM(magi, magj, mags_by_sites, min_sites, total_sites)
    } else if (method == "mutual") {
      res <- mut_infoMM(magi, magj, mags_by_sites, min_sites)
    } else {
      stop("Invalid method")
    }
    
    if (is.null(res)) next
    
    cases_list[[length(cases_list) + 1]] <- res
  }
}
message("\n DONE :)")

# list to data frame
message("\n Preparing output, please wait ...")
cases <- bind_rows(cases_list)

# Correcting multiple testing FDR
if (method == "binomial") {
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
} else if (method == "mutual") {
  write.csv(cases, file = paste0(opt$outdir, 'all_cases.csv'), row.names = FALSE)
}

message("\n Output saved")




