# Network Stats2 (potential edges for ven Diagrams)


#libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
library(optparse)

option_list <- list(
  make_option(c("-m", "--microbial_lineage"), type="character", default="mOTUs_Species_Cluster", help="Name of the microbial lienage"),
  make_option(c("-b", "--bgc_groups"), type="character", default="gcc", help="Name of the grou"),
  make_option(c("-i", "--indir"), type="character", help="Input directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-w", "--workdir"), type="character", help="Working directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
mag_lineage <- opt$microbial_lineage
bgc_group <- opt$bgc_groups

meta_mags <- read.csv(file = paste0(opt$indir, 'metadata.csv'), header = TRUE)
meta_bgcs <- read.csv(file = paste0(opt$indir, 'bgcs_metadata.csv'), header = TRUE)
meta_sites <- read.csv(file = paste0(opt$indir, 'meta_sites.csv'), header = TRUE)

source(paste0(opt$workdir, "functions.R"))

# Count how many microbial lineages and BGC groups are per site
mags_by_sites <- prep_mags(meta_mags, mag_lineage)
bgcs_by_sites<- prep_bgcs(meta_bgcs, bgc_group)

#-------------------------
# MAG - MAG potentials
#-------------------------

# find al de MAGxMAG combiantion that at least co-occur in  1 site

mag_mag_list <- list()
mag_cols <- colnames(mags_by_sites)
mag_cols <- mag_cols[mag_cols != "sites"]
n_mags <- length(mag_cols)
start_time <- Sys.time()

for (i in 1:(n_mags - 1)) {
  
  magi <- mag_cols[i]
  
  if (i %% 100 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    message("\n- Microbial lineages processed:", i)
    message("- Time: ", round(elapsed, 2), " mins ...")
  }
  
  for (j in (i + 1):n_mags) {
    
    magj <- mag_cols[j]
    
    # check if they co-occur 
    comb <- recreate_tableMM(magi, magj, mags_by_sites) # create the table of sites
    shared_sites <- sum(comb[[magi]] > 0 & comb[[magj]] > 0)
    if (shared_sites > 0) {
      mag_mag_list[[length(mag_mag_list) + 1]] <- data.frame(
        MAGi = magi,
        MAGj = magj,
        shared_sites = shared_sites
      )
    }
  }
}
message("\n DONE :)")
# list to data frame
message("\n Preparing MAG-MAG output, please wait ...")
potentials_mm <- bind_rows(mag_mag_list)
write.csv(potentials_mm, file = paste0(opt$outdir, 'potentials_mm.csv'), row.names = FALSE)

#-------------------------
# MAG - BGC potentials
#-------------------------
mag_bgc_list <- list()
#para imprimir avance
counter <- 0
start_time <- Sys.time()

for (col1 in colnames(mags_by_sites)) {
  if (col1 == "sites") next

  counter <- counter + 1
  
  if (counter %% 100 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    message("\n- Microbial lineages processed:", counter)
    message("- Time: ", round(elapsed, 2), " mins ...")
  }
  
  for (col2 in colnames(bgcs_by_sites)) { 
    if (col2 == "sites") next
    
    # check if they co-occur 
    comb <- recreate_tableMB(col1, col2, mags_by_sites, bgcs_by_sites) # create the table of sites
    shared_sites <- sum(comb[[col1]] > 0 & comb[[col2]] > 0)
    if (shared_sites > 0) {
      mag_bgc_list[[length(mag_bgc_list) + 1]] <- data.frame(
        Mags = col1,
        Bgcs = col2,
        shared_sites = shared_sites
      )
    }
  }
}
message("\n DONE :)")

# list to data frame
message("\n Preparing MAG-BGC output, please wait ...")
potentials_mb <- bind_rows(mag_bgc_list)
num_cases <- nrow(potentials_mb)
# filter interactions where the BGC is in the Genome
meta_bgcs <- meta_bgcs %>%
  left_join(meta_mags %>% select(Genome, all_of(mag_lineage)), by = "Genome") # add lineage to bgc table by genome

potentials_mb <- potentials_mb %>% 
  anti_join(meta_bgcs, by = c("Mags" = mag_lineage, "Bgcs" = bgc_group)) 
filt_cases <- num_cases - nrow(potentials_mb)
message("\n NOTE:",filt_cases, " cases where the BGC is in the genome were discarded")
write.csv(potentials_mb, file = paste0(opt$outdir, 'potentials_mb.csv'), row.names = FALSE)
