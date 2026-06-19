##########################################
#### GENERATE PERMUTATIONS ###############
## Andrea Zermeño Díaz ###################
# june-2026 ##############################

# libs
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

# args
option_list <- list(
  make_option(c("-m", "--microbial_lineage"), type="character", default="mOTUs_Species_Cluster", help="Name of the microbial lienage"),
  make_option(c("-b", "--bgc_groups"), type="character", default="gcc", help="Name of the grou"),
  make_option(c("-s", "--minimum_sites"), type="numeric", default=10, help="Minimum number of sites where a group is present"),
  make_option(c("-t", "--temp"), type="character", default="global", help="Range of temperature (max, mid and min)"),
  make_option(c("-w", "--workdir"), type="character", help="Working directory"),
  make_option(c("-i", "--indir"), type="character", help="Input directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-n", "--perms"), type="character", help="Number of permutations"),
)
opt <- parse_args(OptionParser(option_list=option_list))


mag_lineage <- opt$microbial_lineage
bgc_group <- opt$bgc_groups
min_sites <- opt$minimum_sites
temp_r <- opt$temp
method <- opt$method
nperms <- opt$perms

# data load
message("\n Preparing input, please wait ...")
meta_mags <- read.csv(file = paste0(opt$indir, 'metadata.csv'), header = TRUE)
meta_bgcs <- read.csv(file = paste0(opt$indir, 'bgcs_metadata.csv'), header = TRUE)
meta_sites <- read.csv(file = paste0(opt$indir, 'meta_sites.csv'), header = TRUE)
# functions
source(paste0(opt$workdir, "functions.R"))

## Real temperatures #####
# define range
if (temp_r == "low") temp_range <- c(-2, 9)
if (temp_r == "mid") temp_range <- c(10, 20)
if (temp_r == "high") temp_range <- c(21, 35)
if (temp_r != "global") {
  meta_sites <- meta_sites %>%   # filtrar sitios por temperatura
    filter(between(temperature_..C., temp_range[1], temp_range[2]))
}
# filter MAGs y BGCs per sites
meta_mags_real <- meta_mags %>%           
  semi_join(meta_sites, by = "sites")
meta_bgcs_real <- meta_bgcs %>%           
  semi_join(meta_sites, by = "sites")

# count MAGs and 
mags_sites_real <- prep_mags(meta_mags_real, mag_lineage)
bgcs_sites_real <- prep_bgcs(meta_bgcs_real, bgc_group)

# que se 
saveRDS(mags_sites_real, file = file.path(opt$workdir, "perm_0000.rds"))


## Permutated temperatures ########
for (i in seq_len(nperms)) {
  
  sites_perm <- meta_sites
  sites_perm[, -1] <- lapply(sites_perm[, -1, drop = FALSE], sample)
  
  if (temp_r == "low") temp_range <- c(-2, 9)
  if (temp_r == "mid") temp_range <- c(10, 20)
  if (temp_r == "high") temp_range <- c(21, 35)
  if (temp_r != "global") {
    sites_perm <- sites_perm %>%   # filtrar sitios por temperatura
      filter(between(temperature_..C., temp_range[1], temp_range[2]))
  }
  # filter MAGs y BGCs per sites
  meta_mags_perm <- meta_mags %>%           
    semi_join(sites_perm, by = "sites")
  meta_bgcs_perm <- meta_bgcs %>%           
    semi_join(sites_perm, by = "sites")
  
  # count MAGs and 
  mags_sites_perm <- prep_mags(meta_mags_perm, mag_lineage)
  bgcs_sites_perm <- prep_bgcs(meta_bgcs_perm, bgc_group)
  
  
  outfile <- file.path(
    opt$workdir,
    sprintf("perm_%02d.rds", i)
  )
  
  saveRDS(
    mags_perm,
    outfile
  )
}






