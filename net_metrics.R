suppressPackageStartupMessages({
  library(optparse); library(igraph); library(bipartite); library(dplyr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--temp"), make_option("--edges"), make_option("--nodes"),
  make_option("--n_null", type = "integer", default = 100),
  make_option("--n_mod_runs", type = "integer", default = 5),
  make_option("--type_col", default = "type"),      # columna que distingue los dos lados
  make_option("--out"))))

edges <- read.csv(opt$edges, stringsAsFactors = FALSE)
nodes <- read.csv(opt$nodes, stringsAsFactors = FALSE)
if ("node" %in% names(nodes)) nodes <- nodes %>% rename(id = node)
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

tp <- vertex_attr(g, opt$type_col)
V(g)$type <- tp == unique(tp)[1]
stopifnot(is_bipartite(g))
web <- as_biadjacency_matrix(g, types = V(g)$type, sparse = FALSE)
web <- web[rowSums(web) > 0, colSums(web) > 0, drop = FALSE]

nodf <- unname(networklevel(web, index = "NODF"))
Q    <- max(replicate(opt$n_mod_runs, computeModules(web)@likelihood))

nulls     <- nullmodel(web, N = opt$n_null, method = "swap.web")
null_nodf <- sapply(nulls, function(m) unname(networklevel(m, index = "NODF")))
null_Q    <- sapply(nulls, function(m) computeModules(m)@likelihood)

res <- data.frame(temp = opt$temp, n_rows = nrow(web), n_cols = ncol(web),
                  n_links = sum(web > 0),
                  NODF = nodf, NODF_z = (nodf - mean(null_nodf)) / sd(null_nodf),
                  Q = Q,       Q_z    = (Q - mean(null_Q)) / sd(null_Q))
write.csv(res, opt$out, row.names = FALSE)