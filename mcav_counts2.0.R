#TagSeq R script for M.cav FL Keys project
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())
library(arrayQualityMetrics)

#-------------------- Counts ----------------------------
#load counts table
countsmcav <- read.table("mcav_files/mc_feature_counts_buff300_from_saf.tsv", sep = "\t", header = T, row.names = "Geneid")

#remove .fastq.bam from sample names in counts table
countsmcav <- countsmcav %>%
  dplyr::rename_all(
    funs(stringr::str_replace_all(.,".fastq.bam","")) 
  )

#remove non-count columns
countsmcav = select(countsmcav, -c("Chr","Start","End","Length","Strand"))

#load metadata
coldata <- read.table("mcav_files/mcav.metadata.2b_set2", header=F)
names(coldata) <- c('id','age','habitat','admix')

#fix name issue in coldata:
coldata[40,"id"] <- "MA6N"

#change MA1Ob to MA1O in counts
colnames(countsmcav)[colnames(countsmcav) == "MA1Ob"] <- "MA1O"

# Remove hybrid samples and juvenile samples from coldata
coldata_puread <- subset(coldata, age %in% "adult") %>%
  subset(!(admix %in% c("NO1.2","OD2.3","OD2.4")))

#remove individuals in mcavcounts that are not in the coldata list
countsmcav_common <- countsmcav[,colnames(countsmcav) %in% coldata_puread$id]

#remove individuals in coldata that are not in counts_common list
coldata_common <- coldata_puread[coldata_puread$id %in% colnames(countsmcav_common),]

#35 pure lineage, adult mcav samples with 2brad and tagseq data

# Remove genes with cumulative mean read count <5 
min.count <- 5
cum.geneCnt <- apply(countsmcav_common,1,mean)
table(cum.geneCnt >= min.count)
countsmcav_common <- countsmcav_common[cum.geneCnt >= min.count,]
#8546 genes >5 counts

#Run deseq: needed for outlier detection
dds <- DESeqDataSetFromMatrix(countsmcav_common,
                              colData = coldata_common, 
                              design = formula(~habitat + admix))

# Variance stabilizing transformation in order to minimize correlation between the variance and mean
Vsd <- vst(dds)
# creating normalized dataframe
vsd <- assay(Vsd)
# renames the column names
colnames(vsd) <- colnames(countsmcav_common)

# running outlier detection
outliers_adult=ExpressionSet(vsd, AnnotatedDataFrame(as.data.frame(colData(Vsd))))
arrayQualityMetrics(outliers_adult,intgroup=c("habitat","admix"), force=T)

#remove outlier samples: (per quality metrics) MA10D
countsmcav_common <- subset(countsmcav_common, select=-c(MA10D))
coldata_common <- subset(coldata_common, !(id %in% c("MA10D")))

#also removing MA8D (only individual in D3 lineage) - leaves 33 individuals 
countsmcav_common <- subset(countsmcav_common, select=-c(MA8D))
coldata_common <- subset(coldata_common, !(id %in% c("MA8D")))
#habitat: 14 nearshore, 12 offshore, 7 deep
#lineage: 13 N1, 12 O2, 8 D4
#habitat x lineage: 9 nearshore N1, 4 offshore N1, 5 nearshore O2, 7 offshore O2, 1 offshore D3, 7 deep D3

# make the samples the row name for coldata
rownames(coldata_common) <- coldata_common[,1]
coldata_common$id <- NULL

#rename coldata for downstream:
coldata <- coldata_common

#Rerun deseq one last time for remaining 33 samples:
dds <- DESeqDataSetFromMatrix(countsmcav_common,
                              colData = coldata_common, 
                              design = formula(~ habitat + admix))
# Variance stabilizing transformation in order to minimize correlation between the variance and mean
Vsd <- vst(dds)
# creating normalized dataframe
vsd <- assay(Vsd)
# renames the column names
colnames(vsd) <- colnames(countsmcav_common)

#save counts, coldata files
save(countsmcav_common, coldata_common, file = "mcav_files/counts_adults.Rdata")
save(dds, coldata_common, vsd, design, file="mcav_files/vsd_adults.RData")

