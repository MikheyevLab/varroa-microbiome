all_sample_abundance = read.csv("/Users/lazzataibekova/Dropbox (OIST)/Mikheyev Lab/Git/varroa-microbiome/all_sample_abundance.csv", header = T)
library(vegan)
library(ggplot2)

my_data = all_sample_abundance[,1:8786]
row.names(my_data) = all_sample_abundance$species
dissmatrix <- vegdist(my_data, method = "euclidian")
m = monoMDS(dissmatrix, model = "global")
plot(m)
points(m, col = all_sample_abundance$samples, pch = 19)

MDS_xy <- data.frame(m$points)
MDS_xy = data.frame(MDS1 = m$points[,1], MDS2 = m$points[,2])
MDS_xy$samples <- rownames(MDS_xy)
MDS_xy$species <- as.character(all_sample_abundance$species)
MDS_xy$country <- as.character(all_sample_abundance$country)
MDS_xy$host <- as.character(all_sample_abundance$host)

p1 = ggplot(MDS_xy, aes(MDS1, MDS2, color = MDS_xy$host, shape = MDS_xy$species))+
  geom_point()+ 
  coord_equal() + scale_colour_brewer(palette = "Set2") + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) 
p2 = ggplot(MDS_xy, aes(MDS1, MDS2, color = MDS_xy$country, shape = MDS_xy$species))+
  geom_point()+
  coord_equal() + scale_colour_brewer(palette = "Set2") + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) 
ggarrange(p1, p2, legend = "right", ncol = 2, nrow = 1)


data.scores = data.frame(scores(m))
data.scores$Sample = all_sample_abundance[8789]
ggplot(data.scores, aes(MDS1, MDS2, col=Sample))

ggplot() + 
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=species, col=species),alpha=0.5) +  
  coord_equal() +
  theme_bw()

  

(permanova.specieshost <- adonis(dissmatrix ~ species * host, data = binned, permutations = 999, method = "euclidean"))  # ns
(permanova.specieshost <- adonis(dissmatrix ~ species + host, data = binned, permutations = 999, method = "euclidean"))  # ns
(permanova.country <- adonis( , data = binned, permutations = 999, method = "euclidean"))  # sig
(permanova.wol <- adonis(dissmatrix ~ Wolbachia, data = binned, permutations = 999, method = "euclidean"))  # ns