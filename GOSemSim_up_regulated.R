library(GOSemSim)
library(org.Dm.eg.db)
library(data.table)
library(magrittr)
library(tidyverse)
#Reference for Fly
d = godata('org.Dm.eg.db', ont = "BP", computeIC = FALSE)
contrast = read.csv("contrasts_mite_life_stages.csv", header = F)

similarities_upregulated = list()
similarities_downregulated = list()
contrast_names = list()

setwd("GO_enrichment_microbiome_output")
Microbiome_GOs = list.files(pattern = "*.csv") %>% map_df(~read_csv(., col_types = cols(.default = "c")))
GO_Microbiome = as.character(Microbiome_GOs$GO.ID)
setwd("../GO_enrichment_output")
Mite_GOs = list.files(pattern = "*.csv") %>% map_df(~read_csv(., col_types = cols(.default = "c")))
setwd("..")
GO_Mite = as.character(Mite_GOs$GO.ID)

ALL = append(GO_Microbiome, GO_Mite)
GOSEMSIM_matrix = matrix(0, 12, nrow=12)
for (i in 1:(length(contrast$V1)-1))
{
  first_contrast = toString(contrast$V1[i])
  second_contrast = toString(contrast$V2[i])
  filename = paste("./GO_enrichment_output/contrast_", first_contrast, "_", second_contrast, "_upregulated", ".csv", sep = "")
  GO1 = read.csv(filename, header = T)
  GO1ID = as.character(GO1[,1])
  for(j in (i+1):(length(contrast$V1)){
    print(paste(i, j "id"))
    first_contrast = toString(contrast$V1[j])
    second_contrast = toString(contrast$V2[j])
    filename = paste("./GO_enrichment_output/contrast_", first_contrast, "_", second_contrast, "_upregulated", ".csv", sep = "")
    GO2 = read.csv(filename, header = T)
    GO2ID = as.character(GO2[,1])

    GOSEMSIM_matrix[i,j] <- mgoSim(GO1ID, GO2ID, semData = d, measure = "Wang", combine = "BMA")
  }
  print(paste(id, "id"))
  first_contrast = toString(contrast$V1[id])
  second_contrast = toString(contrast$V2[id])
  name = paste(first_contrast, "_", second_contrast,sep = "")
  contrast_names[length(contrast_names)+1] = name

  #UP-REGULATED
  filename = paste("./GO_enrichment_output/contrast_", first_contrast, "_", second_contrast, "_upregulated", ".csv", sep = "")
  GO1 = read.csv(filename, header = T)
  GO1ID = as.character(GO1[,1])

  similarities_upregulated[[length(similarities_upregulated)+1]] = mgoSim(GO1ID, ALL, semData = d, measure = "Wang", combine = "BMA")

  #DOWN-REGULATED
  filename = paste("./GO_enrichment_output/contrast_", first_contrast, "_", second_contrast, "_downregulated", ".csv", sep = "")
  GO1 = read.csv(filename, header = T)
  GO1ID = as.character(GO1[,1])

  similarities_downregulated[[length(similarities_downregulated)+1]] = mgoSim(GO1ID, ALL, semData = d, measure = "Wang", combine = "BMA")
}

names(similarities_downregulated) = as.vector(unlist(contrast_names))
names(similarities_upregulated) = as.vector(unlist(contrast_names))
GOSemSim_table = do.call(rbind, Map(data.frame, A=similarities_upregulated, B=similarities_downregulated))
write.csv(GOSemSim_table, "GOSemSim_table_v2.csv")