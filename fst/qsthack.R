library(vegan)
library(ggplot2)
library(cowplot)
load('mcav_files/vsd_adults.RData')
load("mcav_files/counts_adults.Rdata")
source("https://raw.githubusercontent.com/devanmcg/IntroRangeR/master/11_IntroMultivariate/ordinationsGGplot.R")
library(gridExtra)
library(ggpubr)
library(tidyverse)

#limit vsd to only nearshore and offshore individuals
nvo <- coldata[coldata$habitat != "deep",]

vsd.nvo <- vsd[,rownames(nvo)]

#make dataframe of transformed vsd
datExpr = as.data.frame(t(vsd.nvo))

#calculate capscale distances
dds.pca=capscale(dist(datExpr, method = 'manhattan')~1)
scores.vegan <- as.data.frame(scores(dds.pca,display = "sites", scaling = 1, choices = c(1:10)))

#make a matrix of capscale scores variance stabilized counts using Manhattan distance
scores.vegan.ma <- data.matrix(scores.vegan)

#adonis: analysis of variance of scores.vegan matrix (Manhattan distances of VSD)
perms <- 1000000
set.seed(716863)
adonis2(scores.vegan.ma~habitat,nvo,method = "manhattan",
        permutations = perms)

head(scores.vegan.ma)

#make test sets
head(nvo)
minind = min(table(nvo$habitat))
perc_remove = 1/minind
periter_rem = seq(0, 1, perc_remove)
perms <- 1000000
set.seed(716863)
total_kept = c()
ssq = c()
p_g_F = c()
for(i in 2:(length(periter_rem)-1)){
  near = nvo[sample(rownames(nvo[nvo$habitat == "nearshore",]), size = round((periter_rem[i]*table(nvo$habitat)["nearshore"]))),]
  off = nvo[sample(rownames(nvo[nvo$habitat == "offshore",]), size = round((periter_rem[i]*table(nvo$habitat)["offshore"]))),]
  new_meta = rbind(near, off)
  total_kept = c(total_kept, nrow(near) + nrow(off))
  red.scores = scores.vegan.ma[rownames(scores.vegan.ma) %in% rownames(new_meta),]
  adon_res = adonis2(red.scores~habitat,new_meta,method = "manhattan",
          permutations = perms)
  ssq = c(ssq, adon_res$SumOfSqs[1])
  p_g_F = c(p_g_F, adon_res$`Pr(>F)`[1])
}
results_loop = data.frame(cbind(total_kept, ssq, p_g_F))
results_loop$perc_data = results_loop$total_kept/nrow(nvo)
plot(results_loop$perc_data, results_loop$ssq)
