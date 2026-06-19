################################################
### Permuted temperatures interactions #########
## Andrea Zermeño Díaz ########################
# June-2026 ################################

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
  make_option(c("-p", "--perm"), type="character", help="permutation file")
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
meta_sites <- readRDS(opt$perm) #uno de los sites_perm/real_sites que vienen del primer script
# functions
source(paste0(opt$workdir, "functions.R"))

##### TEMPERATURE ######
# filtrar MAGs y BGCs por sitios
meta_mags <- meta_mags %>%           
  semi_join(meta_sites, by = "sites")
meta_bgcs <- meta_bgcs %>%           
  semi_join(meta_sites, by = "sites")

#### PREP TABLES ####
mags_by_sites <- prep_mags(meta_mags, mag_lineage)
bgcs_by_sites<- prep_bgcs(meta_bgcs, bgc_group)

############### Workflow ##############################

cases_list <- list()
#para imprimir avance
counter <- 0
start_time <- Sys.time()

for (col1 in colnames(mags_by_sites)) {
  if (col1 == "sites") next
  
  # Imprimir progreso
  counter <- counter + 1
  if (counter %% 100 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    message("\n- Microbial lineages processed:", counter)
    message("- Time: ", round(elapsed, 2), " mins ...")
  }
  
  for (col2 in colnames(bgcs_by_sites)) { 
    if (col2 == "sites") next
    
    res <- binomial_MB(col1, col2, mags_by_sites, bgcs_by_sites, min_sites)
    if (is.null(res)) next
    
    cases_list[[length(cases_list) + 1]] <- res
  }
}
message("\n DONE :)")

# list to data frame
message("\n Preparing output, please wait ...")
cases <- bind_rows(cases_list)
num_cases <- nrow(cases)

# filter interactions where the BGC is in the Genome
meta_bgcs <- meta_bgcs %>%
  left_join(meta_mags %>% select(Genome, all_of(mag_lineage)), by = "Genome") # add lineage to bgc table by genome

cases <- cases %>% 
  anti_join(meta_bgcs, by = c("Mags" = mag_lineage, "Bgcs" = bgc_group)) 
filt_cases <- num_cases - nrow(cases)

message("\n NOTE:",filt_cases, " cases where the BGC is in the genome were discarded")

cases <- cases %>%
  mutate(fdr_pval_e = p.adjust(pvalue_e, method = "BH"), 
         fdr_pval_o = p.adjust(pvalue_o, method = "BH"))
occurrence <- cases %>%
  filter(fdr_pval_o <= 0.05)
exclusion <- cases %>%
  filter(fdr_pval_e <= 0.05)

# save the tables 
perm_id <- tools::file_path_sans_ext(basename(opt$perm))
out_all <- file.path(opt$outdir, paste0("all_cases_", perm_id, ".csv"))
out_occ <- file.path(opt$outdir, paste0("oc_filt_", perm_id, ".csv"))
write.csv(cases, out_all, row.names = FALSE)
write.csv(occurrence, out_occ, row.names = FALSE)

message("\n Output saved")