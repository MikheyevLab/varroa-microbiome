library(topGO)
library(sleuth)
library(tidyverse)

#This Function is to build gene2GO map for topGO
readMappings2 = function (file, sep = ",", IDsep = "\t")
{
  b <- read.delim(file = file, header = TRUE, row.names = 1, quote = "",
                  sep = sep, colClasses = "character")
  lst <- by(b[1:length(b)], rownames(b), as.character, simplify = TRUE)
  lst2 = lapply(lst, function(x) x[x != ""])
  return(lst2)
}

#outputs list of genes that have p-value < 0.05
topDiffGenes = function (allScore) {
  return(allScore < 0.05)
}

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

###MAIN###
bacteria <- read_tsv("mondet/gene_calls_nr.bacteria.tsv", col_names = c("contig", "protein_id", "taxname"), col_types = "-c---c-c") %>% filter(grepl("Bacteria", taxname)) %>% mutate(wolbachia = ifelse(grepl("Wolbachia", taxname), "yes", "no")) %>% separate(taxname, sep = ";\\s+", into = c("Bacteria",LETTERS[1:6]), extra = 'drop') %>% mutate(F = sub(";", "", F))
samples <- tibble(id = dir(file.path("mondet/kallisto"))) %>% mutate(path = paste0("mondet/kallisto/", id)) %>%
  left_join(read_tsv("mondet/treatments.tsv", col_types =  cols())) %>% mutate(sample = id, condition = shortName)

contrast = read.csv("contrasts_mite_life_stages.csv", header = F)
gene2GO = readMappings2("GOterm/vdesgoassoc.csv")
t2g = read.csv("IDrna_gene_map.csv", header=T)
for (id in 1:length(contrast$V1))
{
  print(paste(id, "id"))
  first_contrast = toString(contrast$V1[id])
  second_contrast = toString(contrast$V2[id])

  so <- sleuth_prep(dplyr::filter(samples, shortName %in% c(first_contrast, second_contrast)), target_mapping = t2g, aggregation_column = 'ens_gene', gene_mode = T, extra_bootstrap_summary = TRUE)
  so <- sleuth_fit(so, ~condition, 'full')
  so <- sleuth_fit(so, ~1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  cond = colnames(so$fits$full$design_matrix)[2]
  so = sleuth_wt(so, cond)
  results_table = sleuth_results(so, cond)
  #DE_genes = subset(results_table, qval <= 0.05)

  #UP-REGULATED
  up_regulated = subset(results_table, b > 0) #beta postive is up-regulated

  #Preapre inputs for sampleOGData of class topGOdata
  geneList = up_regulated[,2]
  names(geneList) = as.character(up_regulated[,1])
  geneList = sort(geneList, decreasing = TRUE)

  #Build sampleOGdata
  sampleGOdata <- new("topGOdata", description = "Simple session", ontology = "BP", allGenes = geneList, geneSel = topDiffGenes,  nodeSize = 3, annot = annFUN.gene2GO, gene2GO = gene2GO)

  #Fisher's exact test - based on gene counts
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

  #Kolmogorov_smirnov - based on gene scores
  resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

  allGO = usedGO(sampleGOdata)
  allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = length(allGO))

  filename = paste("./GO_enrichment_microbiome_output/contrast_", first_contrast, "_", second_contrast, "_upregulated", ".csv", sep = "")
  write.csv(allRes, file = filename, row.names = F)


  #DOWN-REGULATED
  down_regulated = subset(results_table, b < 0) #beta negative is down-regulated

  #Preapre inputs for sampleOGData of class topGOdata
  geneList = down_regulated[,2]
  names(geneList) = as.character(down_regulated[,1])
  geneList = sort(geneList, decreasing = TRUE)

  #Build sampleOGdata
  sampleGOdata <- new("topGOdata", description = "Simple session", ontology = "BP", allGenes = geneList, geneSel = topDiffGenes,  nodeSize = 3, annot = annFUN.gene2GO, gene2GO = gene2GO)

  #Fisher's exact test - based on gene counts
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

  #Kolmogorov_smirnov - based on gene scores
  resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

  allGO = usedGO(sampleGOdata)
  allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = length(allGO))

  filename = paste("./GO_enrichment_microbiome_output/contrast_", first_contrast, "_", second_contrast, "_downregulated", ".csv", sep = "")
  write.csv(allRes, file = filename, row.names = F)



}
