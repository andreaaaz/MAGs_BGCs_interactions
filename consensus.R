################################################
### Consensus MAG-MAG network ##################
## Andrea Zermeño Díaz #########################
# october-2026 ###############################

# Create a subset "Consensus", which has the shared edges between MAG-MAG y 
# MAG-MAG reconstructed networks, the global shared edges are 11,727 (go to network_QCstats.R)

# libs


# args

# 



# data load
mag_mag <- read.csv("2026-09-interactions08/mOTUs_Species_Cluster/global/oc_filt.csv") # con sitios de alta calidad
mag_mag_r <- read.csv("2026-09-interactions08/mOTUs_Species_Cluster_gcc/global/edges_mm.csv")

# make edges unique identifiers 
make_edge_id <- function(data, node1, node2) { 
  apply( data[, c(node1, node2)], 1, function(x) { 
    paste(sort(x), collapse = "--") } ) 
}

mag_mag$edge <- make_edge_id(mag_mag, "MAGi", "MAGj" )
mag_mag$edge <- make_edge_id(mag_mag, "source", "target")
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
