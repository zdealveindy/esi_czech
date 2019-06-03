# Zeleny D. & Chytry M. (2019): Ecological Specialization Indices for species of the Czech flora. Preslia, 91:93-116. https://doi.org/10.23855/preslia.2019.093
# Supplementary R code 3: Vltava river valley case study
# Author of the R code: David Zeleny (zeleny@ntu.edu.tw, May 2019)


# Required libraries:
library (weimea)  #used version: 0.1.14, devtools::install_github ('zdealveindy/weimea')
library (dplyr)

# Loaded data (from GitHub):
vltava_spe_herbs <- read.delim ('https://raw.githubusercontent.com/zdealveindy/esi_czech/master/vltava_spe_herbs.txt', header = TRUE, check.names = FALSE, row.names = 1)
vltava_nomencl_ESI <- read.delim ('https://raw.githubusercontent.com/zdealveindy/esi_czech/master/vltava_nomenclature_esi.txt', header = TRUE, stringsAsFactors = FALSE)
ESI <- read.delim ('https://raw.githubusercontent.com/zdealveindy/esi_czech/master/Pladias-ESI.txt', header = TRUE, stringsAsFactors = FALSE)
vltava_env <- read.delim ('https://raw.githubusercontent.com/zdealveindy/esi_czech/master/vltava_env.txt', row.names = 1, header = TRUE, stringsAsFactors = FALSE)

vltava_spe_herbs_pa <- vegan::decostand (vltava_spe_herbs, 'pa') # species composition data transformed into presence-absece

# Synchronizing nomenclature between Vltava dataset and the list of ESI values
vltava_ESI_merged <- dplyr::left_join (vltava_nomencl_ESI, ESI[ESI$freq_f >= 10, c('sci.name', 'ESI_f') ], by = c('spec_name_ESI' = 'sci.name'))
vltava_ESI <- vltava_ESI_merged[, 'ESI_f', drop = F]
row.names (vltava_ESI) <- as.character (vltava_ESI_merged$spec_name)

# Calculating fourth corner and community weighted mean analysis between ESI_f and 13 env. variables
# (parallel calculation on the computer with 4 cores - change argument parallel = 4 into lower/higher number if your computer has different number of cores)

cwm_ESI <- cwm (com = vltava_spe_herbs_pa, traits = vltava_ESI$ESI_f)
set.seed (5867)
vltava.cwm <- test_cwm (cwm = cwm_ESI, env = vltava_env, perm = 49999, adjustP = TRUE, p.adjust.method = 'fdr', parallel = 4, test = c('row', 'col', 'max'))
vltava.fco <- test_fourth (cwm = cwm_ESI, env = vltava_env, perm = 49999, adjustP = TRUE, p.adjust.method = 'fdr', parallel = 4, test = 'max')

vltava.cwm.fco.f <- cbind (coef (vltava.cwm)[, c('r', 'P_max', 'P_max_adj')], coef (vltava.fco)[,c('r_ch', 'P_max', 'P_max_adj')])
vltava.cwm.fco.f
write.table (vltava.cwm.fco.f, file = 'clipboard', sep = '\t')



