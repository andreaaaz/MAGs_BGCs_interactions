################################################
### Identifying interactions MAG-BGC ##########
## Andrea Zermeño Díaz #########################
# february-2026 ################################

# libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
library(optparse)


# Args
option_list <- list(
  make_option(c("-m", "--microbial_lineage"), type="character", default="mOTUs_Species_Cluster", help="Name of the microbial lienage"),
  make_option(c("-b", "--bgc_groups"), type="character", default="gcf", help="Name of the grou"),
  make_option(c("-s", "--minimum_sites"), type="numeric", default=10, help="Minimum number of sites where a group is present"),
  make_option(c("-i", "--indir"), type="character", help="Input directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-w", "--workdir"), type="character", help="Working directory"),
  make_option(c("-t", "--temp"), type="character", default="global", help="Range of temperature (max, mid and min)"),
  make_option(c("-e", "--method"), type="character", default="mutual", help="Method to calculate significance")
)
opt <- parse_args(OptionParser(option_list=option_list))

mag_lineage <- opt$microbial_lineage
bgc_group <- opt$bgc_groups
min_sites <- opt$minimum_sites
temp_r <- opt$temp
method <- opt$method

# Run script
# Rscript -m mOTUs_Species_Cluster -b gcf -s 5 -t global -i /mnt/atgc-d3/sur/users/azermeno/exp/MAGs_BGCs_interactions/
# -o /mnt/atgc-d3/sur/users/azermeno/exp/2025-interacions/


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
# filtrar MAGs y BGCa por sitios
meta_mags <- meta_mags %>%           
  semi_join(meta_sites, by = "sites")
meta_bgcs <- meta_bgcs %>%           
  semi_join(meta_sites, by = "sites")

# Count how many microbial lineages and BGC groups are per site
mags_by_sites <- prep_mags(meta_mags, mag_lineage)
bgcs_by_sites<- prep_bgcs(meta_bgcs, bgc_group)


############### Workflow ##############################

# tabla donde se va a guardar informacion de los patrones
cases_list <- list()
#para imprimir avance
counter <- 0
start_time <- Sys.time()

# esto generara 'm x b' tablas de 806 renglones (sitios definidos por station_depth)
# donde 'm' es el numero de microbial lienages
# y 'b' es el numero de grupos de BGCs 


for (col1 in colnames(mags_by_sites)) {
  if (col1 == "sites") next
  
  # Imprimir progreso
  counter <- counter + 1
  if (counter %% 100 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    # percentage <- round((counter / ncol(motus_by_site)) * 100, 2)
    
    message("\n- Microbial lineages processed:", counter)
    message("- Time: ", round(elapsed, 2), " mins ...")
    
  }
  
  for (col2 in colnames(bgcs_by_sites)) { 
    if (col2 == "sites") next
    
    # choosing method
    if (method == "binomial") {
      res <- binomial_MB(col1, col2, mags_by_sites, bgcs_by_sites, min_sites)
    } else if (method == "mutual") {
      res <- mut_infoMB(col1, col2, mags_by_sites, bgcs_by_sites, min_sites)
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
num_cases <- nrow(cases)

# filter interactions where the BGC is in the Genome
meta_bgcs <- meta_bgcs %>%
  left_join(meta_mags %>% select(Genome, all_of(mag_lineage)), by = "Genome") # add lineage to bgc table by genome

cases <- cases %>% 
  anti_join(meta_bgcs, by = c("Mags" = mag_lineage, "Bgcs" = bgc_group)) 
filt_cases <- num_cases - nrow(cases)

message("\n NOTE:",filt_cases, " cases where the BGC is in the genome were discarded")

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




