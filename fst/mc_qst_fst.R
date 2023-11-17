#### calculate gst from gene expression in mcav
library(tidyverse)
library(dplyr)
library(DESeq2)
library(viridis)
library(cowplot)
library(vegan)

#load adult data:
load("mcav_files/vsd_adults.RData")

#first split vsd by sites
near <- coldata[coldata$habitat == "nearshore",]
off <- coldata[coldata$habitat == "offshore",]
deep <- coldata[coldata$habitat == "deep",]

vsd.near <- vsd[,rownames(near)]
vsd.off <- vsd[,rownames(off)]
vsd.deep <- vsd[,rownames(deep)]

#nead to calculate the distance within and between populations

#calculate manhattan distance for each vsd habitat matrix (within pop difference)
vsd.near.man <- as.matrix(dist(vsd.near, method = "manhattan"))
vsd.off.man <- as.matrix(dist(vsd.off, method = "manhattan"))
vsd.deep.man <- as.matrix(vegdist(vsd.deep, method = "manhattan"))

#need to calculate the distance between populations: nvo, nvd, ovd



#

















#average the per-gene manhattan distance for each habitat
qst.near <- rowMeans(vsd.near.man)
qst.off <- rowMeans(vsd.off.man)
qst.deep <- rowMeans(vsd.deep.man)

#load fst data and snp data
load("mcav_files/fst_files/mcavsnps_adults.rdata")

# Pairwise per-locus Fst
fst.files <- list.files('mcav_files/fst_files/outputs_from_tacc/')
fst.ls <- list()
for (i in fst.files) {
  fst.ls[[gsub('\\.fst', '', i)]] <- read.table(paste0('mcav_files/fst_files/outputs_from_tacc/', i))
}

# Count the number of SNPs per gene
snps.per.gene <- group_by(allsnps2, gene) %>% summarise(n_snps = n())

#Fst
fst.log2foldchange.df <- fst.ls[['pno']] %>%
  mutate(V3 = ifelse(V3 < 0, 0, V3)) %>%
  mutate(pos = paste(V1, V2, sep = ':')) %>%
  right_join(snp_coord2, by = 'pos') %>%
  group_by(gene)%>%
  summarise(fst = sum(V3)/sum(V3+V4))
%>%
  left_join(pair.degs, by = 'gene') %>%
  filter(!is.na(log2fold)) %>%
  filter(!is.na(fst))
overallfst <- fst.ls[['po2d3']] %>%
  mutate(V3 = ifelse(V3 < 0, 0, V3)) %>%
  summarise(ofst = sum(V3)/sum(V3+V4)) 

#check if fst is greater/equal/lower on average to deg:
qst.fst <- fst.log2foldchange.df
qst.fst.deg <- as.matrix(qst.fst$log2fold)
#transform log2fold change to be between 0-1
qst.fst$deg <- apply(qst.fst.deg, MARGIN = 2, FUN = function(X) (X-min(X))/diff(range(X)))
#calculate the mean of the difference between fst and deg
qst.fst$difference <- qst.fst$fst - qst.fst$deg
mean(qst.fst$difference)

















#load kinship matrix - row 1 - first individual, row 2 - second individual, K0==J9; K1==J8; K2==J7; rab is the pairwise relatedness
newres <- read.table("mcav_files/newres", header = T)

#load bams
bams <- read_tsv("popgen/mcav_set2/bams_noclones", col_names = 'bam_name') %>% mutate(bam_index = (1:n()) - 1)
bams$bam_name <- gsub("_s2.trim.bt2.bam","",bams$bam_name)
bams$bam_name <- gsub("-r","", bams$bam_name)

#rename index to the bams names
new_res_named = newres %>% as_tibble() %>% 
  left_join(bams %>% dplyr::rename(a = bam_index, bam_a = bam_name)) %>% 
  left_join(bams %>% dplyr::rename(b = bam_index, bam_b = bam_name))

#make kins dataframe
kins <- new_res_named[,c("rab","bam_a","bam_b")]

#load adult data:
load("mcav_files/vsd_adults.RData")

#get list of just adult samples
adults <- rownames(coldata)

#subset to just adult samples
kins2 <- subset(kins, (bam_a %in% adults))
kins2 <- subset(kins2, (bam_b %in% adults))

#create a pairwise matrix
kinmat <- kins2 %>%
  pivot_wider(names_from = c(bam_a), values_from = rab)









#load fst data and snp data
load("mcav_files/fst_files/mcavsnps_adults.rdata")

# Pairwise per-locus Fst
fst.files <- list.files('mcav_files/fst_files/outputs_from_tacc/')
fst.ls <- list()
for (i in fst.files) {
  fst.ls[[gsub('\\.fst', '', i)]] <- read.table(paste0('mcav_files/fst_files/outputs_from_tacc/', i))
}

# Count the number of SNPs per gene
snps.per.gene <- group_by(allsnps2, gene) %>% summarise(n_snps = n())

#plot log2fold change vs. Fst
fst.log2foldchange.df <- fst.ls[['pno']] %>%
  mutate(V3 = ifelse(V3 < 0, 0, V3)) %>%
  mutate(pos = paste(V1, V2, sep = ':')) %>%
  right_join(snp_coord2, by = 'pos') %>%
  group_by(gene)%>%
  summarise(fst = sum(V3)/sum(V3+V4))%>%
  
  left_join(pair.degs, by = 'gene') %>%
  filter(!is.na(log2fold)) %>%
  filter(!is.na(fst))
overallfst <- fst.ls[['po2d3']] %>%
  mutate(V3 = ifelse(V3 < 0, 0, V3)) %>%
  summarise(ofst = sum(V3)/sum(V3+V4)) 

#check if fst is greater/equal/lower on average to deg:
qst.fst <- fst.log2foldchange.df
qst.fst.deg <- as.matrix(qst.fst$log2fold)
#transform log2fold change to be between 0-1
qst.fst$deg <- apply(qst.fst.deg, MARGIN = 2, FUN = function(X) (X-min(X))/diff(range(X)))
#calculate the mean of the difference between fst and deg
qst.fst$difference <- qst.fst$fst - qst.fst$deg
mean(qst.fst$difference)








