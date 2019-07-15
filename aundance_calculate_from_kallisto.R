library(compositions)
library(data.table)
library(dplyr)
setwd("/Users/lazzataibekova/Dropbox (OIST)/Mikheyev Lab/Git/varroa-microbiome/kallisto_geographic")

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                      full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
                    full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}

dir = list.dirs("/Users/lazzataibekova/Dropbox (OIST)/Mikheyev Lab/Git/varroa-microbiome/kallisto_geographic")
filename = file.path(dir, "abundance.tsv")

results = c()
VDnames = c()

VD149 = fread("./VD149/abundance.tsv", header = T, sep = "\t")
VD149_DEG = VD149[target_id %in% DEGs_unique]
#VD149_clr = compositions::clr(VD149$tpm)
VD149_DEG_clr = compositions::clr(VD149_DEG$tpm)
names(VD149_DEG_clr) = VD149_DEG$target_id

my.df =as.data.frame(t(VD149_DEG_clr), col.names = names(VD149_DEG_clr), row.names = VD149)


for (j in seq(along=dir)){
  filename = file.path(dir[[j]], "abundance.tsv")
  VD = fread(filename, header = T, sep = "\t")
  VD_DEG = VD[target_id %in% DEGs_unique]
  #VD149_clr = compositions::clr(VD149$tpm)
  VD_DEG_clr = compositions::clr(VD_DEG$tpm)
  names(VD_DEG_clr) = VD_DEG$target_id

  my.df <- rbind(my.df, VD_DEG_clr)
  VDnames = c(VDnames, dir[[j]])
}

new.df = my.df[-1,]
row.names(new.df) = VDnames


new_info = read.csv("/Users/lazzataibekova/Dropbox (OIST)/Mikheyev Lab/Git/varroa-microbiome/samples.csv", header = T, row.names = 1)
all_sample_abundance = cbind(new.df, new_info)
write.csv(all_sample_abundance, "/Users/lazzataibekova/Dropbox (OIST)/Mikheyev Lab/Git/varroa-microbiome/all_sample_abundance.csv")

