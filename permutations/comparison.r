################################################
### P-values graphs and network statistics #####
## Andrea Zermeño Díaz ########################
# June-2026 ##################################

# libraries
suppressPackageStartupMessages(library(tidyverse))
library(optparse)

# Args
option_list <- list(
  make_option(c("-m", "--microbial_lineage"), type="character", default="mOTUs_Species_Cluster", help="Name of the microbial lienage"),
  make_option(c("-b", "--bgc_groups"), type="character", default="gcc", help="Name of the grou"),
  make_option(c("-s", "--minimum_sites"), type="numeric", default=10, help="Minimum number of sites where a group is present"),
  make_option(c("-t", "--temp"), type="character", default="global", help="Range of temperature (max, mid and min)"),
  make_option(c("-w", "--workdir"), type="character", help="Working directory"),
  make_option(c("-i", "--indir"), type="character", help="Input directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
)
opt <- parse_args(OptionParser(option_list=option_list))
mag_lineage <- opt$microbial_lineage
bgc_group <- opt$bgc_groups
min_sites <- opt$minimum_sites
temp_r <- opt$temp

# recibir todas las permutaciones de cada temperatura que envia nextflow y hacer un csv y guardarlo (por mientras)




