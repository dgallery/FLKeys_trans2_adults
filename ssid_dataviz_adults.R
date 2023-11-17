# Data visualizations with PCA, heatmaps, etc.
#--------------------------------------PCA---------------------------------
library(vegan)
library(ggplot2)
library(cowplot)
load('ssid_files/vsd_adult.RData')
load("ssid_files/counts_adult.Rdata")
source("https://raw.githubusercontent.com/devanmcg/IntroRangeR/master/11_IntroMultivariate/ordinationsGGplot.R")
library(gridExtra)
library(ggpubr)

#make dataframe of transformed vsd
datExpr = as.data.frame(t(vsd))

#calculate capscale distances
dds.pca=capscale(dist(datExpr, method = 'manhattan')~1)
scores.vegan <- as.data.frame(scores(dds.pca,display = "sites", scaling = 1, choices = c(1:4)))
ssid.pca <- as.data.frame(cbind(coldata_common, dds.pca$CA$u))

#create axis labels which include the percent variance of the PCs
choices <- c(1,2)
axis.labels <- ord_labels(dds.pca)[choices]

#####PLOT PCAs by age, site and cluster for the identified PCs of significance
# reorder by site:
ssid.pca$habitat <- factor(ssid.pca$habitat, levels = c("nearshore","offshore","deep"))

#Site PCA: MDS1 & 2
s1_ssid <- ggplot(data = ssid.pca) +
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

s1_ssid

# reorder by admix:
ssid.pca$admix <- factor(ssid.pca$admix, levels = c("S1","D1","D2"))

#Cluster PCA: MDS1 & 2
c1_ssid <- ggplot(data = ssid.pca)+
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab(axis.labels[1]) +
  ylab(axis.labels[2]) +
  geom_point(aes(MDS1, MDS2, fill = admix), shape = 21, size = 5, alpha = 0.8)+
  scale_fill_manual(name = "Ecomorph", labels = c("Shallow 1", "Deep 1", "Deep 2"), 
                    values = c('palegreen1', 'skyblue','dodgerblue4')) +
  labs(fill = "Ecomorph", size = 20) +
  theme(axis.title = element_text(size =16),axis.text = element_text(size = 14),
        legend.title = element_text(size=16), legend.text = element_text(size = 14))

c1_ssid

save(s1_ssid,c1_ssid,file = "ssid_files/ssid_tag_vsd_pca.rdata")

#----------------adonis --------------

#make a matrix of capscale scores variance stabilized counts using Manhattan distance
scores.vegan.ma <- data.matrix(scores.vegan)

#adonis: analysis of variance of scores.vegan matrix (Manhattan distances of VSD)
perms <- 1000000
set.seed(716863)
adonis2(scores.vegan.ma~habitat,coldata_common,method = "manhattan",
        permutations = perms)
#           Df SumOfSqs      R2      F Pr(>F)    
# habitat   2   7566.9 0.45297 19.873  1e-06 ***
# Residual 48   9138.0 0.54703                  
# Total    50  16704.9 1.00000 

set.seed(716863)
adonis2(scores.vegan.ma~admix,coldata_common,method = "manhattan",
        permutations = perms)
#           Df SumOfSqs      R2      F  Pr(>F)  
# admix     2     1665 0.09967 2.6569 0.01192 *
# Residual 48    15040 0.90033                 
# Total    50    16705 1.00000 

set.seed(716863)
adonis2(scores.vegan.ma~habitat*admix,coldata_common,method = "manhattan",
        permutations = perms)
#           Df SumOfSqs      R2       F Pr(>F)    
# habitat   2   7566.9 0.45297 20.1337  1e-06 ***
# admix     2    494.0 0.02957  1.3143 0.2487    
# Residual 46   8644.1 0.51746                   
# Total    50  16704.9 1.00000       

set.seed(716863)
adonis2(scores.vegan.ma~admix*habitat,coldata_common,method = "manhattan",
        permutations = perms)
#           Df SumOfSqs      R2       F   Pr(>F)    
# admix     2   1665.0 0.09967  4.4301 0.000136 ***
# habitat   2   6395.8 0.38287 17.0179    1e-06 ***
# Residual 46   8644.1 0.51746                     
# Total    50  16704.9 1.00000         
