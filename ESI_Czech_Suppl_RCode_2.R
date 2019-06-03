# Zeleny D. & Chytry M. (2019): Ecological Specialization Indices for species of the Czech flora. Preslia, 91:93-116. https://doi.org/10.23855/preslia.2019.093
# Supplementary R code 2: Calculating summary and ploting figures
# Author of the R code: David Zeleny (zeleny@ntu.edu.tw, May 2019)

# install.packages ('ppcor')  # required for calculating partical correlation
# install.packages ('agricolae')  # required for calculating post hoc comparison (HSD.test)

### Load results of calculation, calculate summary and plot figures ----

# Load data from GitHub
load (url ('https://github.com/zdealveindy/esi_czech/blob/master/theta92.r?raw=true')) # load theta92
load (url ('https://github.com/zdealveindy/esi_czech/blob/master/theta73.r?raw=true')) # load theta73
load (url ('https://github.com/zdealveindy/esi_czech/blob/master/theta19.r?raw=true')) # load theta19

# Summary statistics ----
range (theta92$ESI) #2.832964 8.368945
round (quantile (theta92$ESI, prob = c(0, .05, .5, .95, 1)), 2)
#   0%   5%  50%  95% 100% 
# 2.83 3.69 4.97 6.47 8.37

range (theta73$ESI) #2.680678 7.542412
round (quantile (theta73$ESI, prob = c(0, .05, .5, .95, 1)), 2)
#   0%   5%  50%  95% 100% 
# 2.68 3.56 4.93 6.49 7.54 

range (theta19$ESI) #3.154690 7.491409
round (quantile (theta19$ESI, prob = c(0, .05, .5, .95, 1)), 2)
#   0%   5%  50%  95% 100% 
# 3.15 4.08 5.08 6.20 7.49 

range (theta92$occur.freq) # 10 12656
range (theta73$occur.freq) # 10 12442
range (theta19$occur.freq) # 10 5900

# Three the most common species
theta92[order (theta92$occur.freq, decreasing = T), c("sciname_new", "occur.freq")][1:3,]
# sciname_new occur.freq
# 7    Achillea millefolium agg.      12656
# 1493             Urtica dioica      11137
# 1425 Taraxacum sect. Taraxacum      10628

theta73[order (theta73$occur.freq, decreasing = T), c("sciname_new", "occur.freq")][1:3,]
# sciname_new occur.freq
# 7    Achillea millefolium agg.      12442
# 1362 Taraxacum sect. Taraxacum      10043
# 1047        Poa pratensis agg.       9442

theta19[order (theta19$occur.freq, decreasing = T), c("sciname_new", "occur.freq")][1:3,]
# sciname_new occur.freq
# 563 Oxalis acetosella       5900
# 607     Poa nemoralis       5829
# 334   Fagus sylvatica       5352

# number of species for various thresholds
length (theta92$ESI[theta92$occur.freq >= 10]) # 1597
length (theta92$ESI[theta92$occur.freq >= 20]) # 1432
length (theta92$ESI[theta92$occur.freq >= 50]) # 1198
length (theta92$ESI[theta92$occur.freq >= 100])# 1010

length (theta73$ESI[theta73$occur.freq >= 10]) # 1529
length (theta73$ESI[theta73$occur.freq >= 20]) # 1357
length (theta73$ESI[theta73$occur.freq >= 50]) # 1113
length (theta73$ESI[theta73$occur.freq >= 100])# 914

length (theta19$ESI[theta19$occur.freq >= 10]) # 881
length (theta19$ESI[theta19$occur.freq >= 20]) # 747
length (theta19$ESI[theta19$occur.freq >= 50]) # 558
length (theta19$ESI[theta19$occur.freq >= 100])# 415

# Export table with 10 most generalist and 10 most specialist species
theta92.occ <- theta92[theta92$occur.freq >=50,]
genspe92.10 <- theta92.occ[order (theta92.occ$ESI, decreasing = T), c("sciname_new", "ESI", "occur.freq")][c(1:10, (nrow (theta92.occ)-9):nrow (theta92.occ)), ]
write.table (file = 'clipboard', cbind (as.character (genspe92.10$sciname_new), paste (formatC(genspe92.10$ESI,format = 'f', digits = 2), " (", genspe92.10$occur.freq, ")", sep = '')), sep = '\t')

theta73.occ <- theta73[theta73$occur.freq >=50,]
genspe73.10 <- theta73.occ[order (theta73.occ$ESI, decreasing = T), c("sciname_new", "ESI", "occur.freq")][c(1:10, (nrow (theta73.occ)-9):nrow (theta73.occ)), ]
write.table (file = 'clipboard', cbind (as.character (genspe73.10$sciname_new), paste (formatC(genspe73.10$ESI,format = 'f', digits = 2), " (", genspe73.10$occur.freq, ")", sep = '')), sep = '\t')
# 
theta19.occ <- theta19[theta19$occur.freq >=50,]
genspe19.10 <- theta19.occ[order (theta19.occ$ESI, decreasing = T), c("sciname_new", "ESI", "occur.freq")][c(1:10, (nrow (theta19.occ)-9):nrow (theta19.occ)), ]
write.table (file = 'clipboard', cbind (as.character (genspe19.10$sciname_new), paste (formatC(genspe19.10$ESI,format = 'f', digits = 2), " (", genspe19.10$occur.freq, ")", sep = '')), sep = '\t')

# how many species overlaps between forest and non-forest database
sum (theta19$sci.name %in% theta73$sci.name)  # 829 species shared among both datasets

## Plot figures ----
cols <- RColorBrewer::brewer.pal (3, 'PuOr')[c(2,1,3)]

# Fig 1 tiff ----
tiff ('ESI-Fig1.tiff', width = 126, height = 188, units = 'mm', res = 600, pointsize = 12, compression = 'lzw')
par (mfrow = c(3,2), mar = c(4.6, 4.1, 3.1, 2.1))

# Distribution of ESI values ----
hist (theta92$ESI, xlim = c(2,9), xlab = 'Ecological Specialization Index (ESI)', ylab = 'Number of species', main = NA, col = paste (cols[1], '90', sep = ''), breaks = seq (1.5, 9.5, by = .25), las = 1, cex.axis = .9, ylim = c(0,200))
hist (theta73$ESI, add= T, col = paste (cols[2], '90', sep = ''), breaks = seq (1.5, 9.5, by = .25))
hist (theta19$ESI, add = T, col = paste (cols[3], 'CC', sep = ''), breaks = seq (1.5, 9.5, by = .25))
legend ('topright', pch = c(22, 22, 22), col = 1, pt.bg = cols, pt.cex = 2, legend = c('whole', 'non-for.', 'forest'), title = 'Dataset', bty = 'n')
title (main = 'A', adj = 0)

# Relationship of ESI and species frequency ----
plot (ESI ~ occur.freq, theta92, log = 'x', xlab = 'Species frequency', ylab = expression (Ecological~Specialization~Index~(ESI[w])), pch = 16, cex = 0.5, las = 1, col = 'black', cex.axis = .9)
legend ('topright', bty = 'n', legend = bquote (paste (rho==.(formatC (cor (theta92$ESI, theta92$occur.freq, method = 'spearman'), format = 'f', digits = 2)))))
title (main = 'B', adj = 0)

# Relationship of ESI and number of associations per species ----
cor_ESI_Aabs <- cor (theta92$ESI, theta92$nass92.abs, method = 'sp', use = 'pair') #-0.58
pcor_ESI_Aabs_freq <- ppcor::pcor (na.omit (theta92[, c('ESI', 'nass92.abs', 'occur.freq')]), method = 'sp')$estimate['nass92.abs', 'ESI'] # -0.59
plot (ESI ~ nass92.abs, theta92, log = 'x', ylab = expression (Ecological~Specialization~Index~(ESI[w])), xlab = expression (paste ("Number of associations (", A[abs], ")")), pch = 16, cex = 0.5, col = 'black', las = 1, cex.axis = .9)
legend ('topright', bty = 'n', legend = bquote (paste (rho==.(formatC (cor_ESI_Aabs, format = 'f', digits = 2)), ", ", rho[part]==.(formatC (pcor_ESI_Aabs_freq, format = 'f', digits = 2)))))
title (main = 'C', adj = 0)

# correlation of ESI and Aocc, ESI and Arel
cor.test (theta92$ESI, theta92$nass92.occ, method = 'sp') #-0.46
ppcor::pcor (theta92[, c('ESI', 'nass92.occ', 'occur.freq')], method = 'sp')$estimate['nass92.occ', 'ESI'] # -0.59

cor.test (theta92$ESI, theta92$nass92.rel, method = 'sp') #-0.45
ppcor::pcor (na.omit (theta92[, c('ESI', 'nass92.rel', 'occur.freq')]), method = 'sp')$estimate['nass92.rel', 'ESI'] # -0.46


# Specialization and number of habitats (narrowly defined sensu Sadlo et al. 2007) ----
cor_ESI_Hocc <- cor (theta92$ESI, theta92$frequency_narrow_habitats+theta92$frequency_as_optimum_narrow_habitats, use = 'complete', method = 'sp') # -0.47
pcor_ESI_Hocc_freq <- ppcor::pcor (na.omit (cbind (theta92$ESI, theta92$frequency_narrow_habitats+theta92$frequency_as_optimum_narrow_habitats, theta92$occur.freq)), method = 'sp')$estimate[2,1] # -0.43
plot (ESI ~ I(theta92$frequency_narrow_habitats+theta92$frequency_as_optimum_narrow_habitats), theta92, pch = 16, cex = 0.5, col = 'black', las = 1, cex.axis = .9, xlab = expression (paste ('Number of habitats (', H[occ], ')')), ylab = expression (Ecological~Specialization~Index~(ESI[w])))
legend ('topright', bty = 'n', legend = bquote (paste (rho==.(formatC (cor_ESI_Hocc, format = 'f', digits = 2)), ', ', paste (rho[part]==.(formatC (pcor_ESI_Hocc_freq, format = 'f', digits = 2))))))
title (main = 'D', adj = 0)

# correlation with Hopt
cor (theta92$ESI, theta92$frequency_as_optimum_narrow_habitats, use = 'complete', method = 'sp')  # -0.37
ppcor::pcor (na.omit (cbind (theta92$ESI, theta92$frequency_as_otimum_narrow_habitats, theta92$occur.freq)), method = 'sp')$estimate[2,1] # -0.30


Origin2 <- as.character (theta92$origin)
Origin2 <- factor (Origin2, levels = c('orig', 'arch', 'neo'))
ESI2 <- theta92$ESI
bx.inv <- boxplot (ESI2 ~ Origin2, notch = T, outline = F, ylim = c(2,9), names = c('native', 'archaeo', 'neo'), las = 1, xlab = 'Taxon origin', ylab = expression (Ecological~Specialization~Index~(ESI[w])), boxwex = .45, cex.axis = .9)
summary (aov (ESI2 ~ Origin2)) # F =  8.983, P = 0.000132 ***
agricolae:::HSD.test (aov (ESI2 ~ Origin2), trt = "Origin2", group = T, console = T)
text (x = 1:3, y = 8.5, labels = c('a', 'a', 'b'))
text (x = 1:3, y = 2.2, labels = as.vector (table (Origin2)), cex = .7)
title (main = 'E', adj = 0)

# Specialization and conservation status - Czech Republic - not directly used in the paper----
redlist.category <- theta92$red
redlist.category.2 <- substr (as.character (redlist.category), 1,2)
redlist.category.3 <- as.character (redlist.category.2)
redlist.category.3 [redlist.category.3 == ''] <- 'a'
redlist.category.3 <- factor (redlist.category.3, levels = c('A1', 'A2', 'C1', 'C2', 'C3', 'C4', 'a'))
agricolae:::HSD.test (aov (theta92$ESI ~ redlist.category.3), 'redlist.category.3', group = T, console = T)

#IUCN
IUCN <- as.character (theta92$IUCN)
IUCN[IUCN == ''] <- 'LC'
IUCN[IUCN == 'DD'] <- 'NA'
IUCN <- factor (IUCN, levels = c('RE', 'CR', 'EN', 'VU', 'NT', 'LC'))
ESI.w.RE <- theta92$ESI[IUCN != 'RE' | is.na (IUCN)]  # remove RE category (only 4 species)
IUCN.w.RE <- IUCN[IUCN != 'RE' | is.na (IUCN)]
summary (aov (ESI.w.RE ~ IUCN.w.RE)) # F =  88.79, P <2e-16 ***
agricolae:::HSD.test (aov (ESI.w.RE ~ IUCN.w.RE), 'IUCN.w.RE', group = T, console = T)
bx.iucn <- boxplot (ESI.w.RE ~ IUCN.w.RE, notch = T, ylim = c(2,9), names = c('RE', 'CR', 'EN', 'VU', 'NT', 'LC'), xlab = 'Red List (IUCN categories)', ylab = expression (Ecological~Specialization~Index~(ESI[w])), las = 1, outline = F, cex.axis = .88)
points (x = c(1,1,1,1), y = theta92[IUCN == 'RE', 'ESI'][!is.na (theta92[IUCN == 'RE', 'ESI'])], cex = 1.5, pch = 3)
text (x = 2:9, y = 8.5, labels = c('a', 'a', 'a', 'b', 'c'))
text (x = 1:9, y = 2.2, labels = as.vector (table (IUCN)), cex = .7)
title (main = 'F', adj = 0)

dev.off ()

## FIg 2 (Specialization and Ellenberg IV redefined for the Czech Republic) ----

tiff ('ESI-Fig2.tiff', width = 126, height = 188, units = 'mm', res = 600, pointsize = 12, compression = 'lzw')
par (mfrow = c(3,2), mar = c(4.6, 3.1, 3.1, 2.1))
boxplot (ESI ~ light, theta92, las = 1, outline = F, col = c(rep ('white', 8), 'grey'), at = c(2:9, 10.5), names = c(2:9, 'x'))
title (xlab = 'Czech indicator values for light')
title (ylab = expression (Ecological~Specialization~Index~(ESI[w])), line = 1.8)
title (main = 'A', adj = 0)

boxplot (ESI ~ temp, theta92, las = 1, outline = F, col = c(rep ('white', 8), 'grey'), at = c(1:8, 9.5), names = c(1:8, 'x'))
title (xlab = 'Czech indicator values for temperature')
title (ylab = expression (Ecological~Specialization~Index~(ESI[w])), line = 1.8)
title (main = 'B', adj = 0)

boxplot (ESI ~ moist, theta92, las = 1, outline = F, col = c(rep ('white', 12), 'grey'), at = c(1:12, 13.5), names = c(1:12, 'x'))
title (xlab = 'Czech indicator values for moisture')
title (ylab = expression (Ecological~Specialization~Index~(ESI[w])), line = 1.8)
title (main = 'C', adj = 0)

boxplot (ESI ~ react, theta92, las = 1, outline = F, col = c(rep ('white', 9), 'grey'), at = c(1:9, 10.5), names = c(1:9, 'x'))
title (xlab = 'Czech indicator values for soil reaction')
title (ylab = expression (Ecological~Specialization~Index~(ESI[w])), line = 1.8)
title (main = 'D', adj = 0)

boxplot (ESI ~ nutr, theta92, las = 1, outline = F, col = c(rep ('white', 9), 'grey'), at = c(1:9, 10.5), names = c(1:9, 'x'))
title (xlab = 'Czech indicator values for nutrients')
title (ylab = expression (Ecological~Specialization~Index~(ESI[w])), line = 1.8)
title (main = 'E', adj = 0)

boxplot (ESI ~ salin, theta92, las = 1, outline = F)
title (xlab = 'Czech indicator values for salinity')
title (ylab = expression (Ecological~Specialization~Index~(ESI[w])), line = 1.8)
title (main = 'F', adj = 0)
dev.off ()


