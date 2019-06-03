# Zeleny D. & Chytry M. (2019): Ecological Specialization Indices for species of the Czech flora. Preslia, 91:93-116. https://doi.org/10.23855/preslia.2019.093
# Supplementary R code 1: Calculating Ecological Specialization Index
# Author of the R code: David Zeleny (zeleny@ntu.edu.tw, May 2019)

# Note: the script is provided as a reference for the statistical analysis. Since primary data used in this R code
# are not publicly available, this R script is not reproducible. This does not apply to the other two R codes, for which
# data are included at GitHub repository.

# Used libraries:
library (tidyverse)
library (parallel)
# devtools::install_github ('zdealveindy/theta')
library (theta)  # used version 0.7-41

# Custom-built functions:
remove_zerocols <- function (x) x[,colSums (x) > 0, drop = F]  # removes zero columns in species comp data (used in HCR resampling)

no_rel_beta <- function (beta, min_rel, max_rel, beta_max = 1) # predicts number of releves for given betadiversity (sensu Wiser & De Caceres 2013)
{
  beta_min <- (min_rel/max_rel)*beta_max
  slope_k <- max_rel/beta_max
  if (beta >= 0 & beta <= beta_min) no_rel <- min_rel else
    if (beta > beta_min & beta < beta_max) no_rel <- beta*slope_k else
      if (beta >= beta_max & beta <= 1) no_rel <- max_rel
  return (no_rel)
}

remove.blank <- function (x)
{
  unlist (lapply (x, FUN = function (i) 
    {
    splitted <- unlist (strsplit (i, split = ''))
    if (splitted[length (splitted)] == ' ') 
    {while (splitted[length (splitted)] == ' ') splitted <- splitted[-length (splitted)]}
    paste (splitted, collapse = '')
    }))
}

Arel <- function (sitspe, assoc)  # sitspe is species composition data.frame, assoc is vector of associations
{
  assoc <- as.factor (assoc)
  apply (sitspe, 2, FUN = function (foc.spe)
  {
    Ni <- as.vector (table (assoc))
    ni <- as.vector (table (assoc[foc.spe > 0]))
    p_rel_i <- ni/Ni*(1/sum (ni/Ni))  # exponential of Shannon entropy on p_rel_i
    p_rel_i <- p_rel_i[p_rel_i > 0]  # remove zero p_rel_i
    exp (-sum (p_rel_i*log (p_rel_i)))
  })
}

Aabs <- function (sitspe, assoc) # sitspe is species composition data.frame, assoc is vector of associations
{
  assoc <- as.factor (assoc)
  apply (sitspe, 2, FUN = function (foc.spe)
  {
    ni <- as.vector (table (assoc[foc.spe > 0]))
    p_abs_i <- ni*(1/sum (ni))  # exponential of Shannon entropy on p_abs_i
    p_abs_i <- p_abs_i[p_abs_i > 0]  # remove zero p_rel_i
    exp (-sum (p_abs_i*log (p_abs_i)))
  })
}

# Calculate ESI on Czech National Phytosociological Database
# Modify the setwd into the folder where the data are stored:
# setwd ('c:/Users/zeleny/Dropbox/CLANKY/generalists specialists czech vegetation/data/ESI Czech database_v2.0')

# Import database data ----
# whole database
comm92 <- read.table ('92249_Table.txt', sep = '\t', check.names = F, head = T, row.names = 1)
spda92 <- read.table ('92249_SpeciesData.txt', sep = '\t', head = T, fill = T)
head92 <- read.table ('92249_TableHead.txt', sep = '\t', head = T, fill = T)
shhe92 <- read.table ('92249_ShortHead.txt', sep = '\t', head = T, fill = T)

# non-forest database
comm73_0 <- comm92[head92$forest_nonforest == 'nonfor ',]
comm73 <- comm73_0[,colSums (comm73_0) > 0]
spda73 <- spda92[colSums (comm73_0) > 0,]
head73 <- head92[head92$forest_nonforest == 'nonfor ',]
shhe73 <- shhe92[head92$forest_nonforest == 'nonfor ',]

# forest database
comm19_0 <- comm92[head92$forest_nonforest == 'forest ',]
comm19 <- comm19_0[,colSums (comm19_0) > 0]
spda19 <- spda92[colSums (comm19_0) > 0,]
head19 <- head92[head92$forest_nonforest == 'forest ',]
shhe19 <- shhe92[head92$forest_nonforest == 'forest ',]

# ====================================================================================================================
# HCR resampling (each dataset separately) -----
# HCR whole database ----
geog_groups_92 <- unique (shhe92$Group_Number)
no_rel_per_group_92 <- table (shhe92$Group_Number)
comm92_log <- log1p (comm92)

beta_max <- 1
max_rel <- 20
min_rel <- 5

cl <- makeCluster (8)
clusterSetRNGStream(cl, 2308)
clusterExport (cl, varlist = c("comm92_log"))
clusterExport (cl, varlist = c("shhe92", "remove_zerocols", "no_rel_beta", "min_rel", "max_rel", "beta_max"))

res_all_92 <- parLapply (cl, geog_groups_92, fun = function (i) {   #try on i = 981 with 39 plots

  selected_plots <- shhe92$Group_Number == i
  rel <- sum (selected_plots)
  dataset_sel <- remove_zerocols (comm92_log[selected_plots,])
  dataset_sel_dist <- vegan::vegdist (dataset_sel)
  beta <- if (sum (selected_plots) >= 2) mean (dataset_sel_dist) else NA  # calculates mean beta of the plots in the grid
  
  if (rel > min_rel) # do HCR only if the number of plots is larger than min_rel
  {
    no_rel <- round (no_rel_beta (beta = beta, min_rel = min_rel, max_rel = max_rel, beta_max = beta_max))
    if (rel > no_rel) # do HCR only if the number of available plots is larger than the number of required plots
      hcr_selected <- rownames (dataset_sel) [vegclust::hcr (dataset_sel_dist, nout = no_rel, nsampl = 1000)] else
        hcr_selected <- rownames (dataset_sel)
  } else {
    hcr_selected <- rownames (dataset_sel)  # select all rownames in case that min_rel is not reached
  }
  
  beta_hcr <- if (rel > min_rel) mean (vegan::vegdist (remove_zerocols (comm92_log[hcr_selected,]))) else beta
  list (rel = rel, beta = beta, beta_hcr = beta_hcr, hcr_selected = hcr_selected)
})

stopCluster (cl)

hcr_selected_92 <- unlist (sapply (res_all_92, FUN = function (x) x$hcr_selected))
beta_92 <- sapply (res_all_92, FUN = function (x) x$beta)
beta_hcr_92 <- sapply (res_all_92, FUN = function (x) x$beta_hcr)
rel_92 <- sapply (res_all_92, FUN = function (x) x$rel)
rel_hcr_92 <- sapply (res_all_92, FUN = function (x) length (x$hcr_selected))

plot (rel_92 ~ beta_92, ylim = c(0,50))
for (beta in seq (0, 1, length = 1000))
  points (beta, no_rel_beta (beta, min_rel = 5, max_rel = 20, beta_max = 1), pch = '.', cex = 4, col = 'red')
points (rel_hcr_92 ~ beta_hcr_92, col = 'red', pch = 16)

plot (I(beta_hcr_92 - beta_92) ~ I(rel_92-rel_hcr_92), log = 'x')  # after HCR the beta diversity always increases (since the HCR resampled plots are less)

comm92_hcr_0 <- comm92[hcr_selected_92,]
comm92_hcr <- comm92_hcr_0[,colSums (comm92_hcr_0) > 0] 
spda92_hcr <- spda92[colSums (comm92_hcr_0) > 0,]

save ('comm92_hcr', 'spda92_hcr', file = 'comm_and_spda92_hcr.RData' )

# HCR nonforest dataset ----
geog_groups_73 <- unique (shhe73$Group_Number)
no_rel_per_group_73 <- table (shhe73$Group_Number)
comm73_log <- log1p (comm73)

beta_max <- 1
max_rel <- 20
min_rel <- 5

cl <- makeCluster (8)
clusterSetRNGStream(cl, 8254)
clusterExport (cl, varlist = c("comm73_log"))
clusterExport (cl, varlist = c("shhe73", "remove_zerocols", "no_rel_beta", "min_rel", "max_rel", "beta_max"))

res_all_73 <- parLapplyLB (cl, geog_groups_73, fun = function (i) {   
  selected_plots <- shhe73$Group_Number == i
  rel <- sum (selected_plots)
  dataset_sel <- remove_zerocols (comm73_log[selected_plots,])
  dataset_sel_dist <- vegan::vegdist (dataset_sel)
  beta <- if (sum (selected_plots) >= 2) mean (dataset_sel_dist) else NA  # calculates mean beta of the plots in the grid
  
  if (rel > min_rel) # do HCR only if the number of plots is larger than min_rel
  {
    no_rel <- round (no_rel_beta (beta = beta, min_rel = min_rel, max_rel = max_rel, beta_max = beta_max))
    if (rel > no_rel) # do HCR only if the number of available plots is larger than the number of required plots
      hcr_selected <- rownames (dataset_sel) [vegclust::hcr (dataset_sel_dist, nout = no_rel, nsampl = 1000)] else
        hcr_selected <- rownames (dataset_sel)
  } else {
    hcr_selected <- rownames (dataset_sel)  # select all rownames in case that min_rel is not reached
  }
  
  beta_hcr <- if (rel > min_rel) mean (vegan::vegdist (remove_zerocols (comm73_log[hcr_selected,]))) else beta
  list (rel = rel, beta = beta, beta_hcr = beta_hcr, hcr_selected = hcr_selected)
})

stopCluster (cl)

hcr_selected_73 <- unlist (sapply (res_all_73, FUN = function (x) x$hcr_selected))
beta_73 <- sapply (res_all_73, FUN = function (x) x$beta)
beta_hcr_73 <- sapply (res_all_73, FUN = function (x) x$beta_hcr)
rel_73 <- sapply (res_all_73, FUN = function (x) x$rel)
rel_hcr_73 <- sapply (res_all_73, FUN = function (x) length (x$hcr_selected))

plot (rel_73 ~ beta_73, ylim = c(0,50))
for (beta in seq (0, 1, length = 1000))
  points (beta, no_rel_beta (beta, min_rel = 5, max_rel = 20, beta_max = 1), pch = '.', cex = 4, col = 'red')
points (rel_hcr_73 ~ beta_hcr_73, col = 'red', pch = 16)

plot (I(beta_hcr_73 - beta_73) ~ I(rel_73-rel_hcr_73), log = 'x')  # after HCR the beta diversity always increases (since the HCR resampled plots are less)

comm73_hcr_0 <- comm73[hcr_selected_73,]
comm73_hcr <- comm73_hcr_0[,colSums (comm73_hcr_0) > 0]
spda73_hcr <- spda73[colSums (comm73_hcr_0) > 0,]

save ('comm73_hcr', 'spda73_hcr', file = 'comm_and_spda73_hcr.RData' )

# HCR forest dataset ----
geog_groups_19 <- unique (shhe19$Group_Number)
no_rel_per_group_19 <- table (shhe19$Group_Number)
comm19_log <- log1p (comm19)

beta_max <- 1
max_rel <- 20
min_rel <- 5

cl <- makeCluster (8)
clusterSetRNGStream(cl, 7526)
clusterExport (cl, varlist = c("comm19_log"))
clusterExport (cl, varlist = c("shhe19", "remove_zerocols", "no_rel_beta", "min_rel", "max_rel", "beta_max"))

res_all_19 <- parLapplyLB (cl, geog_groups_19, fun = function (i) {   
  selected_plots <- shhe19$Group_Number == i
  rel <- sum (selected_plots)
  dataset_sel <- remove_zerocols (comm19_log[selected_plots,])
  dataset_sel_dist <- vegan::vegdist (dataset_sel)
  beta <- if (sum (selected_plots) >= 2) mean (dataset_sel_dist) else NA  # calculates mean beta of the plots in the grid
  
  if (rel > min_rel) # do HCR only if the number of plots is larger than min_rel
  {
    no_rel <- round (no_rel_beta (beta = beta, min_rel = min_rel, max_rel = max_rel, beta_max = beta_max))
    if (rel > no_rel) # do HCR only if the number of available plots is larger than the number of required plots
      hcr_selected <- rownames (dataset_sel) [vegclust::hcr (dataset_sel_dist, nout = no_rel, nsampl = 1000)] else
        hcr_selected <- rownames (dataset_sel)
  } else {
    hcr_selected <- rownames (dataset_sel)  # select all rownames in case that min_rel is not reached
  }
  
  beta_hcr <- if (rel > min_rel) mean (vegan::vegdist (remove_zerocols (comm19_log[hcr_selected,]))) else beta
  list (rel = rel, beta = beta, beta_hcr = beta_hcr, hcr_selected = hcr_selected)
})

stopCluster (cl)

hcr_selected_19 <- unlist (sapply (res_all_19, FUN = function (x) x$hcr_selected))
beta_19 <- sapply (res_all_19, FUN = function (x) x$beta)
beta_hcr_19 <- sapply (res_all_19, FUN = function (x) x$beta_hcr)
rel_19 <- sapply (res_all_19, FUN = function (x) x$rel)
rel_hcr_19 <- sapply (res_all_19, FUN = function (x) length (x$hcr_selected))

plot (rel_19 ~ beta_19, ylim = c(0,50))
for (beta in seq (0, 1, length = 1000))
  points (beta, no_rel_beta (beta, min_rel = 5, max_rel = 20, beta_max = 1), pch = '.', cex = 4, col = 'red')
points (rel_hcr_19 ~ beta_hcr_19, col = 'red', pch = 16)

plot (I(beta_hcr_19 - beta_19) ~ I(rel_19-rel_hcr_19), log = 'x')  # after HCR the beta diversity always increases (since the HCR resampled plots are less)

comm19_hcr_0 <- comm19[hcr_selected_19,]
comm19_hcr <- comm19_hcr_0[,colSums (comm19_hcr_0) > 0]
spda19_hcr <- spda19[colSums (comm19_hcr_0) > 0,]

save ('comm19_hcr', 'spda19_hcr', file = 'comm_and_spda19_hcr.RData' )

# Calculate Whittaker's beta without outliers for each dataset (parallel computing with 4 cores) ----
theta92.w10.raref.out <- calculate.theta (comm92_hcr, thresh = 10, psample = 10, rarefaction = TRUE, remove.out = TRUE, verbal = T, parallel = TRUE, no.cores = 4)
theta73.w10.raref.out <- calculate.theta (comm73_hcr, thresh = 10, psample = 10, rarefaction = TRUE, remove.out = TRUE, verbal = T, parallel = TRUE, no.cores = 4)
theta19.w10.raref.out <- calculate.theta (comm19_hcr, thresh = 10, psample = 10, rarefaction = TRUE, remove.out = TRUE, verbal = T, parallel = TRUE, no.cores = 4)
# 
# # Save results of calculation for further use ----
save (theta92.w10.raref.out, file = 'theta92.w10.raref.out.r')
save (theta73.w10.raref.out, file = 'theta73.w10.raref.out.r')
save (theta19.w10.raref.out, file = 'theta19.w10.raref.out.r')

# ====================================================================================================================
# Load results of calculation ----
load ('theta92.w10.raref.out.r')
load ('theta73.w10.raref.out.r')
load ('theta19.w10.raref.out.r')

# Use only rarified Whittaker's beta for all calculation:
theta92 <- theta92.w10.raref.out
theta73 <- theta73.w10.raref.out
theta19 <- theta19.w10.raref.out

# Calculate ESI ----
theta92$ESI <- 10-theta92.w10.raref.out$theta
theta73$ESI <- 10-theta73.w10.raref.out$theta
theta19$ESI <- 10-theta19.w10.raref.out$theta


# Rename JUICE names into new nomenclature ---
convert.names <- read.delim ('nomenclature converter_2.0.txt')
# pladias.nomen <- read.delim ('PLADIAS_nomenclature.txt')

theta92 <- merge (theta92, convert.names, by.x = 1, by.y = 1)
theta73 <- merge (theta73, convert.names, by.x = 1, by.y = 1)
theta19 <- merge (theta19, convert.names, by.x = 1, by.y = 1)

# Count number of associations for species ----
# Fix ESY column in head: contains "_unclas" and "_in-mor" values
head92$ESY <- as.character (head92$ESY)
head92$ESY [head92$ESY  %in% c("_unclas", "_in-mor")] <- NA

nass92.occ <- apply (comm92 > 0, 2, FUN = function (x) length (unique (head92$ESY[x][!is.na (head92$ESY[x])])))
nass73.occ <- apply (comm73 > 0, 2, FUN = function (x) length (unique (head73$ESY[x][!is.na (head73$ESY[x])])))
nass19.occ <- apply (comm19 > 0, 2, FUN = function (x) length (unique (head19$ESY[x][!is.na (head19$ESY[x])])))

# A_rel
nass92.rel <- Arel (comm92, head92$ESY)
nass73.rel <- Arel (comm73, head73$ESY)
nass19.rel <- Arel (comm19, head19$ESY)

# A_abs
nass92.abs <- Aabs (comm92, head92$ESY)
nass73.abs <- Aabs (comm73, head73$ESY)
nass19.abs <- Aabs (comm19, head19$ESY)

theta92 <- merge (theta92, data.frame (nass92.occ, nass92.rel, nass92.abs), by.x = 1, by.y = "row.names", all.x = TRUE)
theta73 <- merge (theta73, data.frame (nass73.occ, nass73.rel, nass73.abs), by.x = 1, by.y = "row.names", all.x = TRUE)
theta19 <- merge (theta19, data.frame (nass19.occ, nass19.rel, nass19.abs), by.x = 1, by.y = "row.names", all.x = TRUE)

# Count number of habitats per species ----
habitats <- read.delim ('Sadlo-Chytry-habitats.txt', row.names = 1)
theta92 <- merge(theta92, habitats[,c('frequency_narrow_habitats', 'frequency_as_optimum_narrow_habitats'), drop = F], by.x = 'sciname_new', by.y = 'row.names', all.x = TRUE)
theta73 <- merge(theta73, habitats[,c('frequency_narrow_habitats', 'frequency_as_optimum_narrow_habitats'), drop = F], by.x = 'sciname_new', by.y = 'row.names', all.x = TRUE)
theta19 <- merge(theta19, habitats[,c('frequency_narrow_habitats', 'frequency_as_optimum_narrow_habitats'), drop = F], by.x = 'sciname_new', by.y = 'row.names', all.x = TRUE)


# Ellenberg, redlist and origin in the Czech Rep ----
ell_red_orig <- read.delim ('ell-red-orig.txt', fileEncoding = 'UTF-16LE')

theta92 <- merge(theta92, ell_red_orig, by.x = 'sciname_new', by.y = 'lat_name', all.x = TRUE)
theta73 <- merge(theta73, ell_red_orig, by.x = 'sciname_new', by.y = 'lat_name', all.x = TRUE)
theta19 <- merge(theta19, ell_red_orig, by.x = 'sciname_new', by.y = 'lat_name', all.x = TRUE)

# replace 'x' values from eiv (only for all species in theta92) by 99
CIV_0 <- theta92[, 15:19]
CIV_0 <- apply (CIV_0, 2, as.character)
CIV <- matrix (NA, ncol = ncol (CIV_0), nrow = nrow (CIV_0))
for (co in 1:ncol (CIV_0))
  for (ro in 1:nrow (CIV_0)){
    CIV[ro,co] <- ifelse (any (strsplit (as.character (CIV_0[ro, co]), split = '')[[1]] == 'x'), 99, as.numeric (CIV_0[ro, co]))    
  }

CIV <- sapply (CIV, FUN = function (x) as.numeric (as.character (x)))
theta92[,15:19] <- CIV

# Save data for further calculation ----
save (theta92, file = 'theta92.r')
save (theta73, file = 'theta73.r')
save (theta19, file = 'theta19.r')

# =====================================
# Load data for further calculation ----
load (file = 'theta92.r')
load (file = 'theta73.r')
load (file = 'theta19.r')

# Prepare output table for PLADIAS website ----
merged.all <- theta92[, c("sciname_new", "ESI", "occur.freq")]
names (merged.all) <- c('sci.name', 'ESI_w', 'freq_w')
merged.all <- merge (merged.all, theta73[,c("sciname_new", "ESI", "occur.freq")], by.x = 'sci.name', by.y = 'sciname_new', all.x = TRUE)
names (merged.all)[4:5] <- c('ESI_nf', 'freq_nf')
merged.all <- merge (merged.all, theta19[,c("sciname_new", "ESI", "occur.freq")], by.x = 'sci.name', by.y = 'sciname_new', all.x = TRUE)
names (merged.all)[6:7] <- c('ESI_f', 'freq_f')
write.table (merged.all, file = 'Pladias-ESI.txt', sep = '\t', row.names = F)

