evaluate_pairMM <- function(magi, magj, mags_by_sites, min_sites, total_sites) {
  
  table1 <- mags_by_sites[, c("sites", magi), drop = FALSE]
  table2 <- mags_by_sites[, c("sites", magj), drop = FALSE]
  
  comb <- full_join(table1, table2, by = "sites")
  
  comb[[magi]] <- ifelse(is.na(comb[[magi]]), 0, comb[[magi]])
  comb[[magj]] <- ifelse(is.na(comb[[magj]]), 0, comb[[magj]])
  
  comb <- comb[!(comb[[magi]] == 0 & comb[[magj]] == 0), ]
  
  if (nrow(comb) == 0) return(NULL)
  
  magi_sites <- sum(comb[[magi]] > 0) 
  magj_sites <- sum(comb[[magj]] > 0)
  
  if (magi_sites < min_sites || magj_sites < min_sites) return(NULL)
  
  q <- magi_sites / total_sites
  p <- magj_sites / total_sites
  
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
    pvalue_e = 1 - pbinom(ex_sites, total_sites, pi_e),
    pvalue_o = 1 - pbinom(oc_sites, total_sites, pi_o)
  ))
}

evaluate_pairMB <- function(mag, bgc, m_by_sites, b_by_sites, min_sites, total_sites) {
  
  # create table
  table1 <- m_by_sites[, c("sites", mag), drop = FALSE]
  table2 <- b_by_sites[, c("sites", bgc), drop = FALSE]
  
  temp3 <- full_join(table1, table2, by = "sites")
  
  temp3[[mag]] <- ifelse(is.na(temp3[[mag]]), 0, temp3[[mag]])
  temp3[[bgc]] <- ifelse(is.na(temp3[[bgc]]), 0, temp3[[bgc]])
  
  temp3 <- temp3[!(temp3[[mag]] == 0 & temp3[[bgc]] == 0), ]
  
  if (nrow(temp3) == 0) return(NULL)
  
  mag_sites <- sum(temp3[[mag]] > 0) 
  bgc_sites <- sum(temp3[[bgc]] > 0)
  
  if (mag_sites < min_sites || bgc_sites < min_sites) return(NULL)
  
  q <- mag_sites / total_sites
  p <- bgc_sites / total_sites
  
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
    pvalue_e = 1 - pbinom(ex_sites, total_sites, pi_e),
    pvalue_o = 1 - pbinom(oc_sites, total_sites, pi_o)
  ))
}


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
