library(DESeq2)
library(cowplot)
library(tidyverse)
source('functions/rnaseq_functions.R')

load("ssid_files/counts_adult.Rdata")
load("ssid_files/vsd_adult.RData")

# DESeq Analysis -----------------------------------------------------------

library(BiocParallel)

# The goal of independent filtering is to filter out tests from the procedure that have no,
#or little chance of showing significant evidence, without even looking at their test statistic
# Typically, this results in increased detection power at the same experiment-wide type I error
INDEPENDENT_FILTERING = FALSE

# Running base model for contrast statements 
  #base model = design ~habitat + admix
dds.full <- DESeq(dds, parallel = TRUE)

resultsNames(dds.full)

summary(results(dds.full, independentFiltering = INDEPENDENT_FILTERING))
# out of 6610 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 79, 1.2%
# LFC < 0 (down)     : 89, 1.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)

summary(results(dds.full, contrast = c('habitat', 'nearshore', 'offshore')))
# LFC > 0 (up)       : 1108, 17%
# LFC < 0 (down)     : 1242, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.full, contrast = c('habitat', 'nearshore', 'deep')))
# LFC > 0 (up)       : 1215, 18%
# LFC < 0 (down)     : 1269, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.full, contrast = c('habitat', 'offshore', 'deep')))
# LFC > 0 (up)       : 436, 6.6%
# LFC < 0 (down)     : 413, 6.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

#no S1 adults
summary(results(dds.full, contrast = c('admix', 'S1', 'D1')))
# LFC > 0 (up)       : 79, 1.2%
# LFC < 0 (down)     : 89, 1.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.full, contrast = c('admix', 'S1', 'D2')))
# LFC > 0 (up)       : 1, 0.015%
# LFC < 0 (down)     : 1, 0.015%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.full, contrast = c('admix', 'D1', 'D2')))
# LFC > 0 (up)       : 3, 0.045%
# LFC < 0 (down)     : 3, 0.045%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

# Running DESeq with LRT to compare goodness of fit for a habitat model vs. admix model
# https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
dds.ha <- DESeqDataSetFromMatrix(countssid_common,
                                 colData = coldata_common, 
                                 design = formula(~ habitat + admix))

dds.a=DESeq(dds.ha, test="LRT", reduced = ~habitat, parallel=TRUE) #removing habitat
res_LRT_a <- results(dds.a)
head(res_LRT_a[order(res_LRT_a$padj),])

dds.h=DESeq(dds.ha, test = "LRT", reduced = ~admix, parallel = TRUE) #removing admix
res_LRT_h <- results(dds.h)

#determine which genes differ significantly between admix
sig_res_LRT_a <- res_LRT_a %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  filter(padj < 0.05)
save(sig_res_LRT_a, file = "ssid_files/admix.genes.rdata")

#determine which genes differ significantly between habitat
sig_res_LRT_h <- res_LRT_h %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(padj < 0.05)
save(sig_res_LRT_h, file = "ssid_files/habitat.genes.rdata")

# Creating dataset with design =~habitat (no interaction), to investigate importance of habitat
dds1 <- DESeqDataSetFromMatrix(countssid_common, colData = coldata_common, design = formula(~ habitat))

dds.habitat <- DESeq(dds1, parallel=TRUE)
summary(results(dds.habitat, contrast = c('habitat', 'nearshore', 'offshore')))
# LFC > 0 (up)       : 1099, 17%
# LFC < 0 (down)     : 1231, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.habitat, contrast = c('habitat', 'nearshore', 'deep')))
# LFC > 0 (up)       : 1281, 19%
# LFC < 0 (down)     : 1475, 22%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.habitat, contrast = c('habitat', 'offshore', 'deep')))
# LFC > 0 (up)       : 555, 8.4%
# LFC < 0 (down)     : 539, 8.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

# Creating dataset with design =~admix (no interaction), to investigate importance of admix
dds2 <- DESeqDataSetFromMatrix(countssid_common, colData = coldata_common, design = formula(~ admix))

dds.admix <- DESeq(dds2, parallel=TRUE)
summary(results(dds.admix, contrast = c('admix', 'S1', 'D1')))
# LFC > 0 (up)       : 95, 1.4%
# LFC < 0 (down)     : 133, 2%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.admix, contrast = c('admix', 'S1', 'D2')))
# LFC > 0 (up)       : 1, 0.015%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.admix, contrast = c('admix', 'D1', 'D2')))
# LFC > 0 (up)       : 2, 0.03%
# LFC < 0 (down)     : 2, 0.03%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

#compare differential gene expression of admix+habitat combinations
coldata2 <- coldata_common
coldata2$h.a <- paste0(coldata2$habitat,".",coldata2$admix)

dds3 <- DESeqDataSetFromMatrix(countssid_common, colData = coldata2, design = formula(~ h.a))

dds.combo <- DESeq(dds3, parallel=TRUE)
#compare same habitat, different lineages (only deep has multiple lineages)
summary(results(dds.combo, contrast = c('h.a','deep.S1','deep.D1')))
# LFC > 0 (up)       : 79, 1.2%
# LFC < 0 (down)     : 89, 1.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.combo, contrast = c('h.a','deep.S1','deep.D2')))
# LFC > 0 (up)       : 1, 0.015%
# LFC < 0 (down)     : 1, 0.015%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.combo, contrast = c('h.a','deep.D1','deep.D2'))) 
# LFC > 0 (up)       : 3, 0.045%
# LFC < 0 (down)     : 3, 0.045%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%


#compare same admix, different habitats (only S1 is found in multiple habitats)
summary(results(dds.combo, contrast = c('h.a','nearshore.S1','offshore.S1')))
# LFC > 0 (up)       : 1108, 17%
# LFC < 0 (down)     : 1242, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.combo, contrast = c('h.a','nearshore.S1','deep.S1')))
# LFC > 0 (up)       : 1215, 18%
# LFC < 0 (down)     : 1269, 19%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.combo, contrast = c('h.a','offshore.S1',"deep.S1"))) 
# LFC > 0 (up)       : 436, 6.6%
# LFC < 0 (down)     : 413, 6.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

# saving all models
save(dds.full, dds.h, dds.a, dds.habitat, dds.admix, dds.combo, file="ssid_files/realModels.adults.RData")



#make bar chart of # genes explained by habitat vs lineage
habi <- nrow(sig_res_LRT_h)
lineage <- nrow(sig_res_LRT_a)
num <- c(habi,lineage)
name <- c("Habitat","Lineage")
df <- data.frame(num,name)
ggplot(df, aes(x="",y=num, fill = name))+
  geom_bar(stat = "identity", width=1, color="white")+
  coord_polar("y", start=0) +
  theme_void()+
  guides(fill=guide_legend(title="")) +
  theme(legend.text=element_text(size=16))

