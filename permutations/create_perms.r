##########################################
#### GENERATE PERMUTATIONS ###############
## Andrea Zermeño Díaz ###################
# june-2026 ##############################

# libs
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

# args
option_list <- list(
  make_option(c("-t", "--temp"), type="character", default="global", help="Range of temperature (max, mid and min)"),
  make_option(c("-w", "--workdir"), type="character", help="Working directory"),
  make_option(c("-i", "--indir"), type="character", help="Input directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-n", "--perms"), type="character", help="Number of permutations")
)
opt <- parse_args(OptionParser(option_list=option_list))
temp_r <- opt$temp
nperms <- opt$perms


# data load
message("\n Preparing input, please wait ...")
meta_sites <- read.csv(file = paste0(opt$indir, 'meta_sites.csv'), header = TRUE)


## Real temperatures #####
# define range
if (temp_r == "low") temp_range <- c(-2, 9)
if (temp_r == "mid") temp_range <- c(10, 20)
if (temp_r == "high") temp_range <- c(21, 35)
if (temp_r == "global") {
  real_sites <- meta_sites
} else {
  real_sites <- meta_sites %>%
    filter(between(temperature_..C., temp_range[1], temp_range[2]))
}

real_file <- file.path(opt$outdir, sprintf("real_sites_%s.rds", temp_r))
saveRDS(real_sites, real_file)

## Permuted temperatures ########
for (i in 1:nperms) {
  # permutation of temperatures 
  sites_perm <- meta_sites
  set.seed(i)
  sites_perm$temperature_..C. <- sample(sites_perm$temperature_..C.) 
  # filter sites by temperature 
  if (temp_r == "low") temp_range <- c(-2, 9)
  if (temp_r == "mid") temp_range <- c(10, 20)
  if (temp_r == "high") temp_range <- c(21, 35)
  if (temp_r == "global") {
    sites_perm <- sites_perm
  } else {
    sites_perm <- sites_perm %>%
      filter(between(temperature_..C., temp_range[1], temp_range[2]))
  }
  
  out_file <- file.path(opt$outdir, sprintf("perm_sites_%s_%04d.rds", temp_r, i))
  saveRDS(sites_perm, out_file)
}










