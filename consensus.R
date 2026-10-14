################################################
### Consensus MAG-MAG network ##################
## Andrea Zermeño Díaz #########################
# october-2026 ###############################

# Create a subset "Consensus", which has the shared edges between MAG-MAG y 
# MAG-MAG reconstructed networks, the global shared edges are 11,727 (go to network_QCstats.R)

# libs
suppressPackageStartupMessages(library(tidyverse))

# args
option_list <- list(
  make_option(c("-i", "--indir"), type="character", help="Input directory"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory"),
  make_option(c("-f", "--mm"), type="character", help="file with mm network"),
  make_option(c("-r", "--mmr"), type="character", help="file with mmr network"),
  make_option(c("-w", "--workdir"), type="character", help="Working directory"),
  make_option(c("-t", "--temp"), type="character", default="high", help="Range of temperature (max, mid and min)"))

opt <- parse_args(OptionParser(option_list=option_list))
temp <- opt$temp  

# data load
mag_mag <- read.csv(opt$mm)
mag_mag_r <- read.csv(opt$mmr)
source(paste0(opt$workdir, "functions.R"))

# make edges unique identifiers 
mag_mag$edge <- make_edge_id(mag_mag, "source", "target" )
mag_mag_r$edge <- make_edge_id(mag_mag, "source", "target")
mag_mag_r <- mag_mag_r %>% # ignoring duplicate edges (MAGi <- MAGj, MAGi -> MAGj)
  distinct(edge, .keep_all = TRUE)

# make consensus
consensus <- Reduce(intersect, list(mag_mag$edge, mag_mag_r$edge))
# y recuperamos la tabla
consensus <- data.frame(edge = consensus)
consensus <- consensus %>%
  tidyr::separate(
    edge,
    into = c("MAGi", "MAGj"),
    sep = "--"
  )
# y obtenemos los nodos
nodes_consensus <- data.frame(node = unique(c(consensus$MAGi, consensus$MAGj)))
write.csv(consensus, "mOTUs_Species_Cluster_gcc/global/consensus_edges.csv", row.names = FALSE)
write.csv(nodes_consensus, "mOTUs_Species_Cluster_gcc/global/consensus_nodes.csv", row.names = FALSE)
