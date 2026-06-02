####### Functions ############

#### Binomial test ####

# Make table combined by sites and calculate the p-value of co-exclusion and co-occurrence

### For interactions MAG-MAG ###
binomial_MM <- function(magi, magj, mags_by_sites, min_sites) {
  
  table1 <- mags_by_sites[, c("sites", magi), drop = FALSE]
  table2 <- mags_by_sites[, c("sites", magj), drop = FALSE]
  
  comb <- full_join(table1, table2, by = "sites")
  
  comb[[magi]] <- ifelse(is.na(comb[[magi]]), 0, comb[[magi]])
  comb[[magj]] <- ifelse(is.na(comb[[magj]]), 0, comb[[magj]])
  
  # filtar dobles 0 en la tabla 
  comb <- comb[!(comb[[magi]] == 0 & comb[[magj]] == 0), ]
  n <- nrow(comb)
  # descartar tablas vacias
  if (n == 0) return(NULL)
  # obtener sitios donde aparecen
  magi_sites <- sum(comb[[magi]] > 0) 
  magj_sites <- sum(comb[[magj]] > 0)
  
  # quitar tabla en donde el mag o el bgc este en pocos sitios
  if (magi_sites < min_sites || magj_sites < min_sites) return(NULL)
  
  # foroma menos estricta sitios totales
  # q <- magi_sites / total_sites
  # p <- magj_sites / total_sites
  
  q <- magi_sites/n
  p <- magj_sites/n
  
  
  pi_e <- p - 2*p*q + q
  pi_o <- p * q
  
  ex_sites <- sum((comb[[magi]] == 0 & comb[[magj]] > 0) |
                    (comb[[magi]] > 0 & comb[[magj]] == 0))
  
  oc_sites <- sum((comb[[magi]] > 0) & (comb[[magj]] > 0))
  
  return(list(
    MAGi = magi,
    MAGj = magj,
    MAGi_sites = magi_sites,
    MAGj_sites = magj_sites,
    ex_sites = ex_sites,
    oc_sites = oc_sites,
    q = q,
    p = p,
    pi_exclusion = pi_e,
    pi_occurrence = pi_o,
    pvalue_e = 1 - pbinom(ex_sites, n, pi_e),
    pvalue_o = 1 - pbinom(oc_sites, n, pi_o)
  ))
}

### For MAG-BGC interactions ###
binomial_MB <- function(mag, bgc, m_by_sites, b_by_sites, min_sites) {
  
  table1 <- m_by_sites[, c("sites", mag), drop = FALSE]
  table2 <- b_by_sites[, c("sites", bgc), drop = FALSE]
  
  temp3 <- full_join(table1, table2, by = "sites")
  
  temp3[[mag]] <- ifelse(is.na(temp3[[mag]]), 0, temp3[[mag]])
  temp3[[bgc]] <- ifelse(is.na(temp3[[bgc]]), 0, temp3[[bgc]])
  
  temp3 <- temp3[!(temp3[[mag]] == 0 & temp3[[bgc]] == 0), ]
  
  n <- nrow(temp3)
  if (n == 0) return(NULL)
  
  mag_sites <- sum(temp3[[mag]] > 0) 
  bgc_sites <- sum(temp3[[bgc]] > 0)
  
  if (mag_sites < min_sites || bgc_sites < min_sites) return(NULL)
  
  # q <- mag_sites / total_sites
  # p <- bgc_sites / total_sites
  
  q <- mag_sites/n
  p <- bgc_sites/n
  
  pi_e <- p - 2*p*q + q
  pi_o <- p * q
  
  ex_sites <- sum((temp3[[mag]] == 0 & temp3[[bgc]] > 0) |
                    (temp3[[mag]] > 0 & temp3[[bgc]] == 0))
  
  oc_sites <- sum((temp3[[mag]] > 0) & (temp3[[bgc]] > 0))
  
  return(list(
    Mags = mag, 
    Bgcs = bgc,
    mags_sites = mag_sites, 
    bgcs_sites = bgc_sites,
    ex_sites = ex_sites,
    oc_sites = oc_sites,
    q = q,
    p = p,
    pi_exclusion = pi_e,
    pi_occurrence = pi_o,
    pvalue_e = 1 - pbinom(ex_sites, n, pi_e),
    pvalue_o = 1 - pbinom(oc_sites, n, pi_o)
  ))
}


#### Mutual Information ####

# Create the table and obtain the index of mutual information if the interaction 

### For MAG-BGC interactions ###
mut_infoMB <- function(mag, bgc, m_by_sites, b_by_sites, min_sites){
  
  # crear tabla
  table1 <- m_by_sites[, c("sites", mag), drop = FALSE]
  table2 <- b_by_sites[, c("sites", bgc), drop = FALSE]
  
  temp3 <- full_join(table1, table2, by = "sites")
  
  temp3[[mag]] <- ifelse(is.na(temp3[[mag]]), 0, temp3[[mag]])
  temp3[[bgc]] <- ifelse(is.na(temp3[[bgc]]), 0, temp3[[bgc]])
  
  # saltar las tablas vacias
  if (nrow(temp3) == 0) return(NULL)
  
  # binarizar
  x_bin <- as.integer(temp3[[mag]] > 0)
  y_bin <- as.integer(temp3[[bgc]] > 0)
  
  # probarlo de ambas formas
  if (sum(x_bin) < min_sites || sum(y_bin) < min_sites) return(NULL)
  
  # calcular informacion mutua
  n <- length(x_bin)
  
  x <- x_bin + 1 # indexear
  y <- y_bin + 1
  
  k <- max(c(x, y)) # numero de clasificaciones
  
  Mx <- sparseMatrix(i = 1:n, j = x, x = 1, dims = c(n, k))
  My <- sparseMatrix(i = 1:n, j = y, x = 1, dims = c(n, k))
  
  # Joint
  Pxy <- as.vector(t(Mx) %*% My) / n
  Pxy <- Pxy[Pxy > 0]
  
  Hxy <- -sum(Pxy * log2(Pxy))
  
  # Marginales
  Px <- colMeans(Mx)
  Py <- colMeans(My)
  
  Px <- Px[Px > 0]
  Py <- Py[Py > 0]
  
  Hx <- -sum(Px * log2(Px))
  Hy <- -sum(Py * log2(Py))
  
  # si siempre o nunca aparecen
  if (Hx == 0 || Hy == 0) return(NULL)
  
  MI <- Hx + Hy - Hxy
  # normalizar
  nmi <- sqrt((MI / Hx) * (MI / Hy))
  nmi <- max(0, nmi)
  
  return(tibble(
    Mags = mag,
    Bgcs = bgc,
    NMI = nmi,
    mags_sites = sum(x_bin),
    bgcs_sites = sum(y_bin)
  ))
}

### For MAG-MAG interactions ###
mut_infoMM <- function(magi, magj, m_by_sites, min_sites){
  table1 <- mags_by_sites[, c("sites", magi), drop = FALSE]
  table2 <- mags_by_sites[, c("sites", magj), drop = FALSE]
  
  comb <- full_join(table1, table2, by = "sites")
  
  comb[[magi]] <- ifelse(is.na(comb[[magi]]), 0, comb[[magi]])
  comb[[magj]] <- ifelse(is.na(comb[[magj]]), 0, comb[[magj]])
  
  if (nrow(comb) == 0) return(NULL)
  
  # combertir a 0 y 1
  x_bin <- as.integer(comb[[magi]] > 0)
  y_bin <- as.integer(comb[[magj]] > 0)
  
  if (sum(x_bin) < min_sites || sum(y_bin) < min_sites) return(NULL)
  
  # calcular informacion mutua
  n <- length(x_bin)
  
  x <- x_bin + 1 # indexear
  y <- y_bin + 1
  
  k <- max(c(x, y))
  
  Mx <- sparseMatrix(i = 1:n, j = x, x = 1, dims = c(n, k))
  My <- sparseMatrix(i = 1:n, j = y, x = 1, dims = c(n, k))
  
  # Joint
  Pxy <- as.vector(t(Mx) %*% My) / n
  Pxy <- Pxy[Pxy > 0]
  
  Hxy <- -sum(Pxy * log2(Pxy))
  
  # Marginales
  Px <- colMeans(Mx)
  Py <- colMeans(My)
  
  Px <- Px[Px > 0]
  Py <- Py[Py > 0]
  
  Hx <- -sum(Px * log2(Px))
  Hy <- -sum(Py * log2(Py))
  
  # si siempre o nunca aparecen
  if (Hx == 0 || Hy == 0) return(NULL)
  
  MI <- Hx + Hy - Hxy
  # normalizar
  nmi <- sqrt((MI / Hx) * (MI / Hy))
  nmi <- max(0, nmi)
  
  return(tibble(
    MAGi = magi,
    MAGj = magj,
    NMI = nmi,
    magi_sites = sum(x_bin),
    magj_sites = sum(y_bin)
  ))
}

#### Prepare Tables ####

# Prepare tables for the workflow

prep_mags <- function(meta_mags, mags){
  mags <- sym(mags)
  # data prep
  num_meta <- nrow(meta_mags)
  meta_mags <- meta_mags %>% filter(!is.na(!!mags)) # omit NAs
  fil_meta <- num_meta - nrow(meta_mags)
  message("\n NOTE:", fil_meta, " MAGs were discarded due to missing ", mag_lineage," classification")
  
  # Count the number of MAGs per site and phylogenetic groups
  mags_by_sites <- meta_mags %>%
    count(sites, !!mags) %>%
    pivot_wider(names_from = !!mags, values_from = n, values_fill = 0)
  
  return(mags_by_sites)
}

prep_bgcs <- function(meta_bgcs, bgcs){
  bgcs <- sym(bgcs)
  
  # count the number of BGCs per site and GCF
  bgcs_by_sites <- meta_bgcs %>%
    count(sites, !!bgcs) %>%
    pivot_wider(names_from = !!bgcs, values_from = n, values_fill = 0)
  
  return(bgcs_by_sites)
} 

#### Recreate tables ####
recreate_tableMB <- function(mag, bgc, mags_by_sites, bgcs_by_sites) {
  
  # select the site column and the current column for both tables
  table1 <- mags_by_sites[, c("sites", mag), drop = FALSE]
  table2 <- bgcs_by_sites[, c("
                              sites", bgc), drop = FALSE]
  
  # bind col1 and col2 by site (full_join porque sino se pierden interacciones)
  table_comb <- full_join(table1, table2, by = "sites")
  
  # convertir a 0 los NAs generados por el join
  table_comb[[mag]] <- ifelse(is.na(table_comb[[mag]]), 0, table_comb[[mag]])
  table_comb[[bgc]] <- ifelse(is.na(table_comb[[bgc]]), 0, table_comb[[bgc]])
  
  # filt rows were both are 0
  table_comb <- table_comb[!(table_comb[[mag]] == 0 & table_comb[[bgc]] == 0), ]
  
  return(table_comb)
}

recreate_tableMM <- function(magi, magj, mags_by_sites) {
  
  table1 <- mags_by_sites[, c("sites", magi), drop = FALSE]
  table2 <- mags_by_sites[, c("sites", magj), drop = FALSE]
  
  comb <- full_join(table1, table2, by = "sites")
  
  comb[[magi]] <- ifelse(is.na(comb[[magi]]), 0, comb[[magi]])
  comb[[magj]] <- ifelse(is.na(comb[[magj]]), 0, comb[[magj]])
  
  comb <- comb[!(comb[[magi]] == 0 & comb[[magj]] == 0), ]
  
  return(comb)
}

