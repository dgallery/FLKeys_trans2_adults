#### FST: compare differentially expressed genes with genetic differentiation
library(tidyverse)
library(dplyr)
library(DESeq2)
library(viridis)
library(cowplot)
#########
##### Import datasets #####
load("mcav_files/counts_adults.Rdata")
bams=rownames(coldata)
geno <- read.table('fst/fst_files/sfsSites.geno.gz')
mafs <- read.table('fst/fst_files/sfsSites.mafs.gz', header = T)
genes <- read.table('fst/fst_files/Mcavernosa.maker.coding.3primePlus300.saf', header = T) %>%
  select(gene = GeneID, contig = Chr, start = Start, end = End) %>%
  mutate(contig = as.character(contig), gene = as.character(gene))

#remove leader columns from geno file
cols.rm=ncol(geno) %% length(bams)
snps <- paste(geno[,1], geno[,2], sep = ':')
geno=geno[,-c(1:cols.rm)]

# depending on geno file format, return called genotypes or calculate weighted genotypes based on likelihoods
# split data frame into grouped columns by sample
indCols=ncol(geno)/length(bams)
index=split(1:ncol(geno), rep(1:ncol(geno), each=indCols, length=ncol(geno)))
# calculate weighted genotypes if file matches expected angsd output formats
# -doGeno 6 or 7
if(indCols==2){
  geno=geno[,c(T,F)]
  # -doGeno 8 or 9
} else if(indCols==3){
  geno=data.frame(do.call(cbind, lapply(index, function(i) apply(geno[,i], 1, function(x) sum(x*c(0,1,2))))))
  # -doGeno 10, 11, 12, 13 or 15
} else if(indCols==4 | indCols==5){
  extra=indCols-3
  geno_hard=geno[,seq(1, ncol(geno), by=4)]
  geno=data.frame(do.call(cbind, lapply(index, function(i) apply(geno[,i[-c(1:extra)]], 1, function(x) sum(x*c(0,1,2))))))
} else if(!(indCols %in% c(1:5))){
  stop('Check format of geno file! Incorrect number of columns detected.')
}
names(geno)=bams
tgeno=data.frame(t(geno))
colnames(tgeno) <- snps
# Generate data frame of MAF (-log10 p-value) for functional genes (max MAF if multiple SNPs in gene)
# Identify coordinates of SNP with max MAF within each gene
pb <- txtProgressBar(0,nrow(genes))
Gmax <- c()
snp.maxMAF_coord <- c()

for (i in 1:nrow(genes)) {
  setTxtProgressBar(pb,i)
  gmax <- c()
  snp.maxMAF <- c()
  s <- subset(mafs,chromo==genes$contig[i] & position>=genes$start[i] & position<=genes$end[i])
  if (nrow(s)>0) { #if there are multiple SNPs within a gene...
    # keep the maximum minor allele frequency
    gmax <- c(gmax, max(s$knownEM)) 
    snp.maxMAF <- c(snp.maxMAF, paste(s[base::which.max(s$knownEM),'chromo'], s[base::which.max(s$knownEM),'position'], sep = ':'))
  } else {
    gmax <- c(gmax,0)
    snp.maxMAF <- NA
  }
  Gmax <- rbind(Gmax, gmax)
  snp.maxMAF_coord <- rbind(snp.maxMAF_coord, snp.maxMAF)
}
Gmax <- data.frame(Gmax)
snp_coord <- data.frame(snp.maxMAF_coord) %>% dplyr::rename(pos = 'snp.maxMAF_coord')
Gmax$gene <- genes$gene
snp_coord$gene <- genes$gene
snp_coord2 <- subset(snp_coord, !is.na(pos))

# SNP:gene merge for all SNPs
allsnps <- select(mafs, chromo, position) %>% mutate(chromo = as.character(chromo))
allsnps$gene <- NA
for(i in 1:nrow(allsnps)) {
  gene.match <- which(genes$contig == allsnps$chromo[i] & genes$start <= allsnps$position[i] & genes$end >= allsnps$position[i])
  allsnps$gene[i] <- ifelse(length(gene.match) > 0, genes$gene[gene.match], NA)
}
allsnps2 <- allsnps %>%
  mutate(pos = paste(chromo, position, sep = ':')) %>%
  select(pos, gene) %>%
  filter(!is.na(gene))

#save snps files:
save(snp_coord2,allsnps2,file="mcav_files/fst_files/mcavsnps_adults.rdata")

##### pairwise genes #########
load("Not_in_paper_anymore/fst_files/mcavsnps_adults.rdata")

# Pairwise per-locus Fst
fst.files <- list.files('Not_in_paper_anymore/fst_files/outputs_from_tacc/')
fst.ls <- list()
for (i in fst.files) {
  fst.ls[[gsub('\\.fst', '', i)]] <- read.table(paste0('Not_in_paper_anymore/fst_files/outputs_from_tacc/', i))
}

#deseq models
load("mcav_files/realModels.adults.RData")

#No filtering of DEG or FST by top % fst or adjusted p-val
### Choose pairwise comparison ###
# dds.a = reduced ~habitat (LRT test)
# dds.admix = ~admix
# dds.combo = ~h.a (habitat+admix column)
# dds.full = ~habitat + admix
# dds.h = reduced ~admix (LRT test)
# dds.habitat = ~habitat
pair2compare <- results(dds.a,contrast=c("admix","O2","D3"))

#log2fold change values from DESeq results
deseq.pval <- data.frame(gene = gsub("-RA", "", rownames(pair2compare)), log2fold = pair2compare$log2FoldChange)
pair.degs <- filter(deseq.pval)

# Count the number of SNPs per gene
snps.per.gene <- group_by(allsnps2, gene) %>% summarise(n_snps = n())

#plot log2fold change vs. Fst
fst.log2foldchange.df <- fst.ls[['p23']] %>%
  mutate(V3 = ifelse(V3 < 0, 0, V3)) %>%
  mutate(pos = paste(V1, V2, sep = ':')) %>%
  right_join(snp_coord2, by = 'pos') %>%
  group_by(gene)%>%
  summarise(fst = sum(V3)/sum(V3+V4))%>%
  left_join(pair.degs, by = 'gene') %>%
  filter(!is.na(log2fold)) %>%
  filter(!is.na(fst))
overallfst <- fst.ls[['p23']] %>%
  mutate(V3 = ifelse(V3 < 0, 0, V3)) %>%
  summarise(ofst = sum(V3)/sum(V3+V4)) 

ggplot(fst.log2foldchange.df, aes(y = log2fold, x = fst))+
  ylab(bquote(log[2]~ foldchange))+
  xlab(bquote(F[ST]))+
  geom_hex(bins = 30)+
  scale_fill_viridis(trans="log10") +
  #ggtitle("Nearshore v. Offshore Ecomorph") +
  theme_cowplot()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = overallfst$ofst)



########fst by genes in modules#####
topfst <- function(x){
  goinput = read.csv(paste0("mcav_go/",x,".csv"))
  top1 = inner_join(goinput,fst.log2foldchange.df)
  top2 = top1[order(top1$Fish_kME, decreasing = TRUE),]
  return(top2[1:100,])
}
 
dg_fst_top = topfst("darkgrey")
dt_fst_top = topfst("darkturquoise")
lc_fst_top = topfst("lightcyan")
lg_fst_top = topfst("lightgreen")
gy_fst_top = topfst("greenyellow")
mb_fst_top = topfst("midnightblue")


plot(density(mb_fst_top$fst), col = "midnightblue")
lines(density(lg_fst_top$fst), col = "lightgreen")
lines(density(lc_fst_top$fst), col = "lightcyan")
lines(density(dg_fst_top$fst), col = "darkgrey")
lines(density(dt_fst_top$fst), col = "darkturquoise")
lines(density(gy_fst_top$fst), col="greenyellow")


######qq plots ######

# function to draw ngenes at random ntimes, order, compute mean and sd 
# (returns a list of 2 items: mean and sd)
meanq=function(x, ntimes, ngenes=100){
  dd=list()
  for (i in 1:ntimes) {
    dx=sample(x,ngenes)
    dx=dx[order(dx)] 
    dd[[i]]=dx
  }
  dd=do.call(cbind,dd)
  return(list(apply(dd,1,mean),apply(dd,1,sd)))
}

# assuming you have a vector of all your Fsts called Fst:
#make vector of all fsts
Fst = fst.log2foldchange.df$fst

# logit-transform all your Fsts first:
logit <- function(x,add=0.01){
  x[x<0]=0
  x=x+add
  return(log(x/(1-x)))
}

Fst = logit(Fst)

# these are Fsts for top 100 genes in each module
# (replace with actual code to extract those, this is just faking data to make sure code works)
toplg = lg_fst_top$fst
toplc = lc_fst_top$fst
topdg = dg_fst_top$fst
topdt = dt_fst_top$fst
topgy = gy_fst_top$fst
topmb = mb_fst_top$fst

#log transform top 100 genes in modules
toplg2 = logit(toplg)
toplc2 = logit(toplc)
topdg2 = logit(topdg)
topdt2 = logit(topdt)
topgy2 = logit(topgy)
topmb2 = logit(topmb)


# sampling 100 random genes 100 times, 
# computing mean and SD of Fst distribution in such samples
nullq=meanq(Fst,ntimes=100,ngenes=100)

# reordering module Fsts 
toplg2 = toplg2[order(toplg2)]
toplc2 = toplc2[order(toplc2)]
topdg2 = topdg2[order(topdg2)]
topdt2 = topdt2[order(topdt2)]
topgy2 = topgy2[order(topgy2)]
topmb2 = topmb2[order(topmb2)]

# putting together data frame for ggplotting
f2q=data.frame(cbind(null=nullq[[1]],toplg2,toplc2,topdg2,topdt2,topgy2,topmb2))
f2q.s=data.frame(cbind(mean0=nullq[[1]]),sd0=nullq[[2]],stack(f2q))
names(f2q.s)[3:4]=c("Fst","comparison")
f2q.s$sd0[f2q.s$comparison!="null"]=0

# qqplot with ribbon showing 1.96*SD of null
ggplot(f2q.s,aes(mean0,Fst,group=comparison))+
  geom_ribbon(aes(ymax=mean0+1.96*sd0,ymin=mean0-1.96*sd0),fill="grey90")+
  geom_line(aes(color=comparison))+
  theme_void()+
  scale_color_manual(values = c("red","lightgreen","lightcyan","darkgrey",
                                "darkturquoise","greenyellow","midnightblue"), 
                     labels = c("Null","Lightgreen","Lightcyan","Darkgrey",
                                "Darkturquoise","Greenyellow","Midnightblue"),
                     name = "Comparison")

# We hope to see lines for modules getting above the grey ribbon in the top-right corner

#density plot of f2q
library(reshape)
f2q.melt <- melt(f2q)
ggplot(f2q.melt, aes(x = value, fill = variable))+
  geom_density(alpha=0.25)+
  scale_fill_manual(values = c("red","lightgreen","lightcyan","darkgrey","darkturquoise","greenyellow","midnightblue"))
  
#########linear regression between fst and deg:#######
cor.test(fst.log2foldchange.df$fst, fst.log2foldchange.df$log2fold)

####### determine how many and which genes are outliers for both FST and DEG ### 

#find highest 95% of fst values
ecdf_fst <- ecdf(fst.log2foldchange.df$fst)
topfst <- which(round(ecdf_fst(fst.log2foldchange.df$fst),3) == 0.90)
fst_thresh <- min(fst.log2foldchange.df$fst[topfst])
plot(ecdf(fst.log2foldchange.df$fst))
abline(v=fst_thresh)
abline(h=0.9)

#find the highest 97.5% and lowest 2.5% DEG 
ecdf_log <- ecdf(fst.log2foldchange.df$log2fold)
toplog <- which(round(ecdf_log(fst.log2foldchange.df$log2fold),3) == 0.95)
botlog <- which(round(ecdf_log(fst.log2foldchange.df$log2fold), 3) == .05)

highlog_thresh <- min(fst.log2foldchange.df$log2fold[toplog])
lowlog_thresh <- max(fst.log2foldchange.df$log2fold[botlog])

plot(ecdf(fst.log2foldchange.df$log2fold))
abline(v=highlog_thresh)
abline(v=lowlog_thresh)

abline(h=0.95)
abline(h = 0.05)

#ID outliers
outliers <- fst.log2foldchange.df %>%
  filter(fst >= fst_thresh) %>%
  filter(log2fold >= highlog_thresh | log2fold <= lowlog_thresh)
nrow(outliers)

#check if fst is greater/equal/lower on average to deg:
qst.fst <- fst.log2foldchange.df

qst.fst.deg <- as.matrix(qst.fst$log2fold)

#transform log2fold change to be between 0-1
qst.fst$deg <- apply(qst.fst.deg, MARGIN = 2, FUN = function(X) (X-min(X))/diff(range(X)))

#calculate the mean of the difference between fst and deg
qst.fst$difference <- qst.fst$fst - qst.fst$deg
mean(qst.fst$difference)

#nearshore vs offshore: -0.446256 fst is slightly higher than qst
#nearshore vs deep: -0.5263465 fst is lower than qst
#offshore vs deep: -0.5178728 fst is lower than qst
