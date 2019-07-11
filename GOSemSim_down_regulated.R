library(GOSemSim)
library(org.Dm.eg.db)
library(vegan)

#Reference for Fly
d = godata('org.Dm.eg.db', ont = "BP", computeIC = FALSE)
contrast = read.csv("contrasts_mite_life_stages.csv", header = F)

#MITE
GOSEMSIM_matrix_mite = matrix(0, 12, nrow=12)

#DOWN-REGULATED
for(i in 1:(length(contrast$V1)-1))
{
  first_contrast = toString(contrast$V1[i])
  second_contrast = toString(contrast$V2[i])
  filename = paste("./GO_enrichment_output/contrast_", first_contrast, "_", second_contrast, "_downregulated", ".csv", sep = "")
  GO1 = read.csv(filename, header = T)
  GO1ID = as.character(GO1[,1])

  for(j in (i+1):length(contrast$V1))
  {
    print(paste(i, j))
    first_contrast = toString(contrast$V1[j])
    second_contrast = toString(contrast$V2[j])
    filename = paste("./GO_enrichment_output/contrast_", first_contrast, "_", second_contrast, "_downregulated", ".csv", sep = "")
    GO2 = read.csv(filename, header = T)
    GO2ID = as.character(GO2[,1])
    GOSEMSIM_matrix_mite[i,j] <- mgoSim(GO1ID, GO2ID, semData = d, measure = "Wang", combine = "BMA")
  }
  name = paste(first_contrast, "_", second_contrast,sep = "")
}

#MICROBIOME
GOSEMSIM_matrix_microbiome = matrix(0, 12, nrow=12)

#DOWN-REGULATED
for(i in 1:(length(contrast$V1)-1))
{
  first_contrast = toString(contrast$V1[i])
  second_contrast = toString(contrast$V2[i])
  filename = paste("./GO_enrichment_microbiome_output/contrast_", first_contrast, "_", second_contrast, "_downregulated", ".csv", sep = "")
  GO1 = read.csv(filename, header = T)
  GO1ID = as.character(GO1[,1])

  for(j in (i+1):(length(contrast$V1)))
  {
    print(paste(i, j))
    first_contrast = toString(contrast$V1[j])
    second_contrast = toString(contrast$V2[j])
    filename = paste("./GO_enrichment_microbiome_output/contrast_", first_contrast, "_", second_contrast, "_downregulated", ".csv", sep = "")
    GO2 = read.csv(filename, header = T)
    GO2ID = as.character(GO2[,1])
    GOSEMSIM_matrix_microbiome[i,j] <- mgoSim(GO1ID, GO2ID, semData = d, measure = "Wang", combine = "BMA")
  }
  name = paste(first_contrast, "_", second_contrast,sep = "")
}
microb_dist = as.dist(t(GOSEMSIM_matrix_microbiome))
mite_dist = as.dist(t(GOSEMSIM_matrix_mite))
mantel(mite_dist, microb_dist)

write.csv(as.data.frame(GOSEMSIM_matrix_microbiome), "GOSEMSIM_matrix_microbiome_down.csv")
write.csv(as.data.frame(GOSEMSIM_matrix_mite), "GOSEMSIM_matrix_mite_down.csv")
