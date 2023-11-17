library(DESeq2)
library(cowplot)
library(tidyverse)
source('functions/rnaseq_functions.R')

load("mcav_files/counts_adults.Rdata")
load("mcav_files/vsd_adults.RData")

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
# out of 8546 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.012%
# outliers [1]       : 4, 0.047%
# low counts [2]     : 0, 0%
# (mean count < 0)

summary(results(dds.full, contrast = c('habitat', 'nearshore', 'offshore')))
# LFC > 0 (up)       : 756, 8.8%
# LFC < 0 (down)     : 786, 9.2%
# outliers [1]       : 4, 0.047%
# low counts [2]     : 0, 0%
# (mean count < 5)

summary(results(dds.full, contrast = c('habitat', 'nearshore', 'deep')))
# LFC > 0 (up)       : 11, 0.13%
# LFC < 0 (down)     : 9, 0.11%
# outliers [1]       : 4, 0.047%
# low counts [2]     : 6625, 78%
# (mean count < 79)

summary(results(dds.full, contrast = c('habitat', 'offshore', 'deep')))
# LFC > 0 (up)       : 13, 0.15%
# LFC < 0 (down)     : 4, 0.047%
# outliers [1]       : 4, 0.047%
# low counts [2]     : 6460, 76%
# (mean count < 32)

summary(results(dds.full, contrast = c('admix', 'N1', 'O2')))
# LFC > 0 (up)       : 17, 0.2%
# LFC < 0 (down)     : 23, 0.27%
# outliers [1]       : 4, 0.047%
# low counts [2]     : 0, 0%
# (mean count < 1)

summary(results(dds.full, contrast = c('admix', 'N1', 'D3')))
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.012%
# outliers [1]       : 3, 0.035%
# low counts [2]     : 0, 0%
# (mean count < 1)

summary(results(dds.full, contrast = c('admix', 'O2', 'D3')))
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.012%
# outliers [1]       : 4, 0.047%
# low counts [2]     : 0, 0%
# (mean count < 1)


# Running DESeq with LRT to compare goodness of fit for a habitat model vs. admix model
# https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
dds.ha <- DESeqDataSetFromMatrix(countsmcav_common,
                                 colData = coldata, 
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
save(sig_res_LRT_a, file = "mcav_files/admix.genes.rdata")

#determine which genes differ significantly between habitat
sig_res_LRT_h <- res_LRT_h %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(padj < 0.05)
save(sig_res_LRT_h, file = "mcav_files/habitat.genes.rdata")

# Creating dataset with design =~habitat (no interaction), to investigate importance of habitat
dds1 <- DESeqDataSetFromMatrix(countsmcav_common, colData = coldata, design = formula(~ habitat))

dds.habitat <- DESeq(dds1, parallel=TRUE)
summary(results(dds.habitat, contrast = c('habitat', 'nearshore', 'offshore')))
# LFC > 0 (up)       : 926, 11%
# LFC < 0 (down)     : 989, 12%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.habitat, contrast = c('habitat', 'nearshore', 'deep')))
# LFC > 0 (up)       : 1172, 14%
# LFC < 0 (down)     : 1066, 12%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.habitat, contrast = c('habitat', 'offshore', 'deep')))
# LFC > 0 (up)       : 586, 6.9%
# LFC < 0 (down)     : 423, 4.9%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

# Creating dataset with design =~admix (no interaction), to investigate importance of admix
dds2 <- DESeqDataSetFromMatrix(countsmcav_common, colData = coldata, design = formula(~ admix))

dds.admix <- DESeq(dds2, parallel=TRUE)
summary(results(dds.admix, contrast = c('admix', 'N1', 'O2')))
# LFC > 0 (up)       : 28, 0.33%
# LFC < 0 (down)     : 35, 0.41%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.admix, contrast = c('admix', 'N1', 'D3')))
# LFC > 0 (up)       : 800, 9.4%
# LFC < 0 (down)     : 622, 7.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

summary(results(dds.admix, contrast = c('admix', 'O2', 'D3')))
# LFC > 0 (up)       : 378, 4.4%
# LFC < 0 (down)     : 295, 3.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

#compare differential gene expression of admix+habitat combinations
coldata2 <- coldata
coldata2$h.a <- paste0(coldata2$habitat,".",coldata2$admix)

dds3 <- DESeqDataSetFromMatrix(countsmcav_common, colData = coldata2, design = formula(~ h.a))

dds.combo <- DESeq(dds3, parallel=TRUE)
#compare same habitat, different lineages
summary(results(dds.combo, contrast = c('h.a','nearshore.N1','nearshore.O2')))
# LFC > 0 (up)       : 1, 0.012%
# LFC < 0 (down)     : 8, 0.094%
# outliers [1]       : 1, 0.012%
# low counts [2]     : 0, 0%

summary(results(dds.combo, contrast = c('h.a','offshore.N1','offshore.O2')))
# LFC > 0 (up)       : 9, 0.11%
# LFC < 0 (down)     : 8, 0.094%
# outliers [1]       : 1, 0.012%
# low counts [2]     : 1989, 23%

summary(results(dds.combo, contrast = c('h.a','offshore.N1','offshore.D3'))) #only 1 offshore D3 individual
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 1, 0.012%
# low counts [2]     : 0, 0%

summary(results(dds.combo, contrast = c('h.a','offshore.O2','offshore.D3'))) #only 1 offshore D3 individual
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 3, 0.035%
# outliers [1]       : 1, 0.012%
# low counts [2]     : 0, 0%

#deep only has 1 lineage

#compare same admix, different habitats
summary(results(dds.combo, contrast = c('h.a','nearshore.N1','offshore.N1')))
# LFC > 0 (up)       : 179, 2.1%
# LFC < 0 (down)     : 280, 3.3%
# outliers [1]       : 1, 0.012%
# low counts [2]     : 1823, 21%

summary(results(dds.combo, contrast = c('h.a','nearshore.O2','offshore.O2')))
# LFC > 0 (up)       : 330, 3.9%
# LFC < 0 (down)     : 214, 2.5%
# outliers [1]       : 1, 0.012%
# low counts [2]     : 995, 12%

summary(results(dds.combo, contrast = c('h.a','offshore.D3',"deep.D3"))) #only 1 offshore D3 individual
# LFC > 0 (up)       : 12, 0.14%
# LFC < 0 (down)     : 3, 0.035%
# outliers [1]       : 1, 0.012%
# low counts [2]     : 6958, 81%

# saving all models
save(dds.full, dds.h, dds.a, dds.habitat, dds.admix, dds.combo, file="mcav_files/realModels.adults.RData")

#create go input for results of genes from habitat and admix:
load("mcav_files/admix.genes.rdata")
load("mcav_files/habitat.genes.rdata")

admix_genes <- sig_res_LRT_a[,c(1,3)]
habi_genes <- sig_res_LRT_h[,c(1,3)]

write.csv(admix_genes, file = "mcav_go/admix_genes.csv",row.names=F,quote=F)
write.csv(habi_genes, file = "mcav_go/habi_genes.csv",row.names=F,quote=F)


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
 

