library(GOSemSim)
library(org.Dm.eg.db)
#Reference for Fly
d = godata('org.Dm.eg.db', ont = "BP", computeIC = FALSE)
contrast = read.csv("contrasts_mite_life_stages.csv", header = F)
similarities_upregulated = c()
similarities_downregulated = c()

for (id in 1:length(contrast$V1))
{
  first_contrast = toString(contrast$V1[id])
  second_contrast = toString(contrast$V2[id])

  #UP-REGULATED
  filename = paste("./GO_enrichment_microbiome_output/contrast_", first_contrast, "_", second_contrast, "_upregulated", ".csv", sep = "")
  GO1 = read.csv(filename, header = T)
  GO1ID = as.character(GO1[,1])

  filename = paste("./GO_enrichment_output/contrast_", first_contrast, "_", second_contrast, "_upregulated", ".csv", sep = "")
  GO2 = read.csv(filename, header = T)
  GO2ID = as.character(GO2[,1])

  similarities_upregulated[id] = mgoSim(GO1ID, GO2ID, semData = d, measure = "Wang", combine = "BMA")

  #DOWN-REGULATED
  filename = paste("./GO_enrichment_microbiome_output/contrast_", first_contrast, "_", second_contrast, "_downregulated", ".csv", sep = "")
  GO1 = read.csv(filename, header = T)
  GO1ID = as.character(GO1[,1])

  filename = paste("./GO_enrichment_output/contrast_", first_contrast, "_", second_contrast, "_downregulated", ".csv", sep = "")
  GO2 = read.csv(filename, header = T)
  GO2ID = as.character(GO2[,1])

  similarities_downregulated[id] = mgoSim(GO1ID, GO2ID, semData = d, measure = "Wang", combine = "BMA")
}