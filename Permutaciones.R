################################################
### Permutations for Mutual Information MAG-BGC ##########
## Andrea Zermeño Díaz #########################
# april-2026 ################################

# para saber si el  NMI es por azar o si si esta detectando algo
# se haran permutaciones en los sitios de la tabla mags_by_sites unicamente

# libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
library(parallel)
library(optparse)
library(ggplot2)

# Args
option_list <- list(
  make_option(c("-m", "--microbial_lineage"), type="character", default="mOTUs_Species_Cluster", help="Name of the microbial lienage"),
  make_option(c("-b", "--bgc_groups"), type="character", default="gcc", help="Name of the grou"),
  make_option(c("-s", "--minimum_sites"), type="numeric", default=10, help="Minimum number of sites where a group is present"),
  make_option(c("-t", "--temp"), type="character", default="low", help="Range of temperature (max, mid and min)"),
  make_option(c("-i", "--indir"), type="character", help="Input directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory")
  )
opt <- parse_args(OptionParser(option_list=option_list))

mag_lineage <- opt$microbial_lineage
bgc_group <- opt$bgc_groups
min_sites <- opt$minimum_sites
temp_r <- opt$temp


#### DATA LOAD ####
message("\n Preparing input, please wait ...")

meta_mags <- read.csv(file = paste0(opt$indir, 'metadata.csv'), header = TRUE)
meta_bgcs <- read.csv(file = paste0(opt$indir, 'bgcs_metadata.csv'), header = TRUE)
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
# filtrar MAGs y BGCs por sitios
meta_mags <- meta_mags %>%           
  semi_join(meta_sites, by = "sites")
meta_bgcs <- meta_bgcs %>%           
  semi_join(meta_sites, by = "sites")


#### PREP TABLES ####
mags_by_sites <- prep_mags(meta_mags, mag_lineage)
bgcs_by_sites<- prep_bgcs(meta_bgcs, bgc_group)

############### Workflow ##############################

#### Permutaciones ####

# funcion de permutaciones
run_one_perm <- function(i, mags_by_sites, bgcs_by_sites, min_sites) {
  # set up
  start_perm <- Sys.time()
  counter <- 0
  message("\nRunning permutation ", i)
  # para el resultado real
  if (i == 1) {
    mags_perm <- mags_by_sites
  } else {
  # para las permutaciones
    mags_perm <- mags_by_sites
    mags_perm[, -1] <-
      mags_perm[sample(nrow(mags_perm)), -1]
  }
  
  cases_list <- list()
  
  for (col1 in colnames(mags_perm)) {
    if (col1 == "sites") next
    
    # progreso interno
    counter <- counter + 1
    if (counter %% 100 == 0) {
      elapsed <- difftime(Sys.time(), start_perm, units = "mins")
      message("\n- Columns processed: ", counter, " | Time: ", round(elapsed, 2), " mins")
    }
    
    for (col2 in colnames(bgcs_by_sites)) {
      if (col2 == "sites") next
      
      res <- mut_infoMB(col1, col2, mags_perm, bgcs_by_sites, min_sites)
      
      if (!is.null(res)) {
        cases_list[[length(cases_list) + 1]] <- res
      }
    }
  }
  
  bind_rows(cases_list) %>%
    mutate(perm = i, type = ifelse(i == 1, "real", "perm"))
}

### Parallelization ###

message("--------PERMUTATIONS-------")
# set up 
ncores <- 11
n_perm <- 11
perm_results <- mclapply(
  1:n_perm,
  run_one_perm,
  mags_by_sites = mags_by_sites,
  bgcs_by_sites = bgcs_by_sites,
  min_sites = min_sites,
  mc.cores = ncores
)

message("Saving data, please wait...")
cases_perm <- bind_rows(perm_results)
message("Done :)")


#### Plot distributions ####



#### Save info ####

write.csv(cases_perm, paste0(opt$outdir, "cases_perm.csv"), row.names = FALSE)
