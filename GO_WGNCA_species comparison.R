#wgcna species comparison

#### renaming data to keep track of species #######
#load mcav data:
load("mcav_files/networkdata_signed.adult.RData")
load("mcav_files/wgcnaData.adult.RData")
load("mcav_files/wgcna_adults_input.Rdata")

#rename all data:
mc_datExpr <- datExpr
mc_datt <- datt
mc_datTraits <- datTraits
mc_genetree <- geneTree
mc_MEs <- MEs
mc_traits <- traits
mc_moduleColors <- moduleColors
mc_moduleLabels <- moduleLabels

#load ssid data:
load("ssid_files/networkdata_signed.adult.ssid.RData")
load("ssid_files/wgcnaData.adult.RData")
load("ssid_files/wgcna_adults_input.Rdata")

#rename ssid data:
ss_datExpr <- datExpr
ss_datt <- datt
ss_datTraits <- datTraits
ss_genetree <- geneTree
ss_MEs <- MEs
ss_traits <- traits
ss_moduleColors <- moduleColors
ss_moduleLabels <- moduleLabels

#save species labeled data:
save(mc_datExpr,mc_datt,mc_datTraits,mc_genetree,mc_MEs,mc_traits,mc_moduleColors,mc_moduleLabels,
     ss_datExpr,ss_datt,ss_datTraits,ss_genetree,ss_MEs,ss_traits,ss_moduleColors,ss_moduleLabels,
     file = "species_wgcna.rdata")

##### running comparisons ####
library(WGCNA)
library(flashClust)
library(ape)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()
library(DESeq2)
library(ggplot2)
library(vegan)
library(tidyverse)
library(dplyr)

load("species_wgcna.rdata")

#first comparison: mc light cyan and ss royal blue

#lightcyan:
mcav <- merge(mc_MEs, mc_traits, by = "row.names") 
mcav_lc <- mcav[,-c(3:7)]
mcav_lc$species = "M. cavernosa"
mcav_lc$ME = mcav_lc$MElightcyan

#royalblue
ssid <- merge(ss_MEs,ss_traits, by = "row.names")
ssid_rb <- ssid[,-c(2:4,6:7)]
ssid_rb$species = "S.siderea"
ssid_rb$ME = ssid_rb$MEroyalblue

#combine to one dataframe:
mclc_ssrb <- merge(mcav_lc,ssid_rb, all = T)

#order by habitat
mclc_ssrb$habitat = factor(mclc_ssrb$habitat, levels = c("nearshore","offshore","deep"))

ggplot(data=mclc_ssrb, aes(x = habitat, y = ME, color = species))+
  geom_boxplot(width = 0.75, outlier.shape = NA)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.title = element_text(size = 16))+
  scale_color_discrete(name = "Species", 
                     labels = c(expression(italic("M. cavernosa")), expression(italic("S. siderea"))))+
  scale_y_continuous(name = expression("Module Eigengene"))+
  scale_x_discrete(name = expression("Habitat"), labels = c("Nearshore","Offshore","Deep"))+
  theme(legend.text.align = 0, )+
  ggtitle(bquote('Lightcyan and Royalblue,'~ R^2~'= 0.24, p < 0.05'))
                     

#second comparison: mc greenyellow and ss green

#greenyellow:
mcav_gy <- mcav[,-c(2:5,7)]
mcav_gy$species = "M. cavernosa"
mcav_gy$ME = mcav_gy$MEgreenyellow

#green
ssid_g <- ssid[,-c(2,4:7)]
ssid_g$species = "S.siderea"
ssid_g$ME = ssid_g$MEgreen

#combine to one dataframe:
mcgy_ssg <- merge(mcav_gy,ssid_g, all = T)

#order by habitat
mcgy_ssg$habitat = factor(mclc_ssrb$habitat, levels = c("nearshore","offshore","deep"))

ggplot(data=mcgy_ssg, aes(x = habitat, y = ME, color = species))+
  geom_boxplot(width = 0.75, outlier.shape = NA)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.title = element_text(size = 16))+
  scale_color_discrete(name = "Species", 
                       labels = c(expression(italic("M. cavernosa")), expression(italic("S. siderea"))))+
  scale_y_continuous(name = expression("Module Eigengene"))+
  scale_x_discrete(name = expression("Habitat"), labels = c("Nearshore","Offshore","Deep"))+
  theme(legend.text.align = 0)+
  ggtitle(bquote('Greenyellow and Green,'~ R^2~'= 0.27, p < 0.05'))



#third comparison: mc greenyellow and ss lightyellow

#greenyellow:
mcav_gy <- mcav[,-c(2:5,7)]
mcav_gy$species = "M. cavernosa"
mcav_gy$ME = mcav_gy$MEgreenyellow

#green
ssid_ly <- ssid[,-c(2:3,5:7)]
ssid_ly$species = "S.siderea"
ssid_ly$ME = ssid_ly$MElightyellow

#combine to one dataframe:
mcgy_ssly <- merge(mcav_gy,ssid_ly, all = T)

#order by habitat
mcgy_ssly$habitat = factor(mcgy_ssly$habitat, levels = c("nearshore","offshore","deep"))

ggplot(data=mcgy_ssly, aes(x = habitat, y = ME, color = species))+
  geom_boxplot(width = 0.75, outlier.shape = NA)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        plot.title = element_text(size = 16))+
  scale_color_discrete(name = "Species", 
                       labels = c(expression(italic("M. cavernosa")), expression(italic("S. siderea"))))+
  scale_y_continuous(name = expression("Module Eigengene"))+
  scale_x_discrete(name = expression("Habitat"), labels = c("Nearshore","Offshore","Deep"))+
  theme(legend.text.align = 0)+
  ggtitle(bquote('Greenyellow and Lightyellow,'~ R^2~'= 0.20, p < 0.05'))
