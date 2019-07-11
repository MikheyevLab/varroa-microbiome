library(GOSemSim)
library(org.Dm.eg.db)
#Reference for Fly
d = godata('org.Dm.eg.db', ont = "BP", computeIC = FALSE)
contrast = read.csv("contrasts_mite_life_stages.csv", header = F)

similarities_upregulated = list()
similarities_downregulated = list()
contrast_names = list()

for (id in 1:length(contrast$V1))
{
  print(paste(id, "id"))
  first_contrast = toString(contrast$V1[id])
  second_contrast = toString(contrast$V2[id])

  #UP-REGULATED
  filename = paste("./GO_enrichment_microbiome_output/contrast_", first_contrast, "_", second_contrast, "_upregulated", ".csv", sep = "")
  GO1 = read.csv(filename, header = T)
  GO1ID = as.character(GO1[,1])

  filename = paste("./GO_enrichment_output/contrast_", first_contrast, "_", second_contrast, "_upregulated", ".csv", sep = "")
  GO2 = read.csv(filename, header = T)
  GO2ID = as.character(GO2[,1])

  name = paste(first_contrast, "_", second_contrast,sep = "")
  contrast_names[length(contrast_names)+1] = name
  similarities_upregulated[[length(similarities_upregulated)+1]] = mgoSim(GO1ID, GO2ID, semData = d, measure = "Wang", combine = "BMA")

  #DOWN-REGULATED
  filename = paste("./GO_enrichment_microbiome_output/contrast_", first_contrast, "_", second_contrast, "_downregulated", ".csv", sep = "")
  GO1 = read.csv(filename, header = T)
  GO1ID = as.character(GO1[,1])

  filename = paste("./GO_enrichment_output/contrast_", first_contrast, "_", second_contrast, "_downregulated", ".csv", sep = "")
  GO2 = read.csv(filename, header = T)
  GO2ID = as.character(GO2[,1])

  similarities_downregulated[[length(similarities_downregulated)+1]] = mgoSim(GO1ID, GO2ID, semData = d, measure = "Wang", combine = "BMA")
}

names(similarities_downregulated) = as.vector(unlist(contrast_names))
names(similarities_upregulated) = as.vector(unlist(contrast_names))
GOSemSim_table = do.call(rbind, Map(data.frame, A=similarities_upregulated, B=similarities_downregulated))
write.csv(GOSemSim_table, "GOSemSim_table.csv")