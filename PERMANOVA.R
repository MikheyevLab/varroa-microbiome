library(readxl)
library(vegan)
#library(MASS)
library(compositions)

# Relative abundance by groups
#binned <- read_excel("/Users/lazzataibekova/Dropbox (OIST)/varroa microbiome/stephanie results/Taxonomy_all.xlsx", sheet = "relative abundance")
binned = read.csv("/Users/lazzataibekova/Dropbox (OIST)/Mikheyev Lab/Git/varroa-microbiome/all_sample_abundance.csv", header = T)
binned$samples <- as.factor(binned$samples)
binned$species <- as.factor(binned$species)
binned$host <- as.factor(binned$host)
binned$country <- as.factor(binned$country)
binned$Wolbachia <- as.factor(binned$Wolbachia)
str(binned)



# PERMANOVA
(permanova.specieshost <- adonis(dissmatrix ~ species * host, data = binned, permutations = 999, method = "euclidean"))  # ns
(permanova.specieshost <- adonis(dissmatrix ~ species + host, data = binned, permutations = 999, method = "euclidean"))  # ns
(permanova.country <- adonis( , data = binned, permutations = 999, method = "euclidean"))  # sig
(permanova.wol <- adonis(dissmatrix ~ Wolbachia, data = binned, permutations = 999, method = "euclidean"))  # ns

