# 2bRAD maps -------------------------------
library(ggplot2)
library(maptools)

if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "gshhg-bin-2.3.7/gshhs_f.b"
sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(-83.5, -80), ylim = c(24.2, 26.7)) %>%
  fortify()

#inset map
inset_map <- ggplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_polygon(data = sf1, aes(x=long, y = lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  coord_fixed(xlim = c(-82.5, -80), ylim = c(24.2, 25.7), expand = 0)+
  geom_rect(aes(xmin=-81.6,xmax=-81.42,ymin=24.48,ymax=24.62),fill = NA, color= "red")
inset_map

#zoom map
sample_meta <- as.data.frame(matrix(NA, nrow = 6, ncol = 5))
colnames(sample_meta) <- c("Region", "Reef", "Latitude", "Longitude","Transect")
sample_meta$Region <- "Florida"
sample_meta$Reef <- c("Jaap Reef", "Maryland Shoals Rockpile", "Maryland Shoals Hump","Carly's Patch","Hasluns' Head","Looe Key Deep")
sample_meta$Latitude <- c(24.58570,24.5217166666667,24.49659, 24.60712,24.55320,24.54098)
sample_meta$Longitude <- c(-81.5826166666667,-81.57675,-81.56746,-81.42942,-81.43759,-81.41413)
sample_meta$Transect <- c(2,2,2,1,1,1)

library(leaflet)
library(htmltools)
url = "https://upload.wikimedia.org/wikipedia/commons/c/c0/Location_dot_black.svg"
dot = makeIcon(url,url,24,24)
leaflet(sample_meta) %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  setView( -81.5, 24.545, zoom = 11) %>%
  addMarkers(~Longitude,~Latitude, popup = ~as.character(Reef), icon = dot, label = ~htmlEscape(Reef))


# Data visualizations with PCA, heatmaps, etc.
#--------------------------------------PCA---------------------------------
library(vegan)
library(ggplot2)
library(cowplot)
load('mcav_files/vsd_adults.RData')
load("mcav_files/counts_adults.Rdata")
source("https://raw.githubusercontent.com/devanmcg/IntroRangeR/master/11_IntroMultivariate/ordinationsGGplot.R")
library(gridExtra)
library(ggpubr)

#make dataframe of transformed vsd
datExpr = as.data.frame(t(vsd))

#calculate capscale distances
dds.pca=capscale(dist(datExpr, method = 'manhattan')~1)
scores.vegan <- as.data.frame(scores(dds.pca,display = "sites", scaling = 1, choices = c(1:4)))
mcav.pca <- as.data.frame(cbind(coldata, dds.pca$CA$u))

#create axis labels which include the percent variance of the PCs
choices <- c(1,2)
axis.labels <- ord_labels(dds.pca)[choices]

#####PLOT PCAs by age, site and cluster for the identified PCs of significance
# reorder by site:
mcav.pca$habitat <- factor(mcav.pca$habitat, levels = c("nearshore","offshore","deep"))

#Site PCA: MDS1 & 2
s1_mcav <- ggplot(data = mcav.pca) +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab(axis.labels[1]) +
  ylab(axis.labels[2]) +
  geom_point(aes(MDS1, MDS2, fill = habitat), shape = 21, size =5, alpha = 0.8) +
  scale_fill_manual(name = "Habitat", labels = c("Nearshore", "Offshore", "Deep"), 
                     values = c('#004b23','#8b0a50','#002a52')) +
  labs(fill = "Habitat", size = 20) +
  theme(axis.title = element_text(size =16),axis.text = element_text(size = 14),
        legend.title = element_text(size=16), legend.text = element_text(size = 14))

s1_mcav

# reorder by admix:
mcav.pca$admix <- factor(mcav.pca$admix, levels = c("N1","O2","D3"))

#Cluster PCA: MDS1 & 2
c1_mcav <- ggplot(data = mcav.pca)+
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab(axis.labels[1]) +
  ylab(axis.labels[2]) +
  geom_point(aes(MDS1, MDS2, fill = admix), shape = 21, size = 5, alpha = 0.8)+
  scale_fill_manual(name = "Ecomorph", labels = c("Nearshore", "Offshore", "Deep 1", "Deep 2"), 
                     values = c('palegreen1','plum1','skyblue','dodgerblue4')) +
  labs(fill = "Ecomorph", size = 20) +
  theme(axis.title = element_text(size =16),axis.text = element_text(size = 14),
        legend.title = element_text(size=16), legend.text = element_text(size = 14))

c1_mcav

save(s1_mcav,c1_mcav,file = "mcav_files/mcav_tag_vsd_pca.rdata")

#----------------adonis --------------

#make a matrix of capscale scores variance stabilized counts using Manhattan distance
scores.vegan.ma <- data.matrix(scores.vegan)

#adonis: analysis of variance of scores.vegan matrix (Manhattan distances of VSD)
perms <- 1000000
set.seed(716863)
adonis2(scores.vegan.ma~habitat,coldata,method = "manhattan",
       permutations = perms)
#           Df SumOfSqs      R2      F Pr(>F)    
# habitat   2   8782.4 0.50291 15.175  1e-06 ***
# Residual 30   8680.9 0.49709                  
# Total    32  17463.3 1.00000    

set.seed(716863)
adonis2(scores.vegan.ma~admix,coldata,method = "manhattan",
       permutations = perms)
#          Df SumOfSqs      R2     F Pr(>F)    
# admix     2   4828.8 0.27651 5.733  2e-06 ***
# Residual 30  12634.4 0.72349                 
# Total    32  17463.3 1.00000  

set.seed(716863)
adonis2(scores.vegan.ma~habitat*admix,coldata,method = "manhattan",
       permutations = perms)
#Terms added sequentially (first to last)
#               Df SumOfSqs      R2       F Pr(>F)    
# habitat        2   8782.4 0.50291 15.3232  1e-06 ***
# admix          2    624.0 0.03573  1.0887 0.3727    
# habitat:admix  1    319.4 0.01829  1.1147 0.3630    
# Residual      27   7737.4 0.44307                   
# Total         32  17463.3 1.00000            

set.seed(716863)
adonis2(scores.vegan.ma~admix*habitat,coldata,method = "manhattan",
       permutations = perms)
#               Df SumOfSqs      R2      F Pr(>F)    
# admix          2   4828.8 0.27651 8.4252  1e-06 ***
# habitat        2   4577.5 0.26212 7.9867  1e-06 ***
# admix:habitat  1    319.4 0.01829 1.1147  0.363    
# Residual      27   7737.4 0.44307                  
# Total         32  17463.3 1.00000           
