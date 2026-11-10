suppressPackageStartupMessages({ library(optparse); library(igraph); library(dplyr) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--t1"), make_option("--e1"), make_option("--n1"),
  make_option("--t2"), make_option("--e2"), make_option("--n2"),
  make_option("--out"))))

load_net <- function(ef, nf) {
  edges <- read.csv(ef, stringsAsFactors = FALSE)
  nodes <- read.csv(nf, stringsAsFactors = FALSE)
  if ("node" %in% names(nodes)) nodes <- nodes %>% rename(id = node)
  graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
}

edge_keys <- function(g) {
  el <- as_edgelist(g, names = TRUE)
  data.frame(a = pmin(el[,1], el[,2]), b = pmax(el[,1], el[,2]), stringsAsFactors = FALSE)
}
key <- function(k) paste(k$a, k$b, sep = "||")

sorensen <- function(x, y) {
  a <- length(intersect(x, y)); b <- length(setdiff(x, y)); c <- length(setdiff(y, x))
  if ((2*a + b + c) == 0) return(NA_real_)
  (b + c) / (2*a + b + c)
}

g1 <- load_net(opt$e1, opt$n1); g2 <- load_net(opt$e2, opt$n2)
shared <- intersect(V(g1)$name, V(g2)$name)
k1 <- edge_keys(g1); k2 <- edge_keys(g2)
only_shared <- function(k) k[k$a %in% shared & k$b %in% shared, ]

WN <- sorensen(key(k1), key(k2))
OS <- sorensen(key(only_shared(k1)), key(only_shared(k2)))
S  <- sorensen(V(g1)$name, V(g2)$name)

res <- data.frame(t1 = opt$t1, t2 = opt$t2, S = S, OS = OS, ST = WN - OS, WN = WN,
                  prop_rewiring = OS / WN)
write.csv(res, opt$out, row.names = FALSE)