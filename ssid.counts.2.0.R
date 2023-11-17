#TagSeq R script for S.Sid FL Keys project
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())
library(arrayQualityMetrics)

#-------------------- Counts ----------------------------
#load counts table
countssid <- read.table("ssid_files/allcounts_ssid.txt", sep = "\t", header = T)

#set gene names column to row names and remove non-count columns
rownames(countssid) <- countssid[,1]
countssid[,1:2] <- NULL

#load metadata
coldata <- read.table("ssid_files/ssid.metadata.2b_set2", header=F)
names(coldata) <- c('id','age','habitat','admix')

#fix name issues in coldata:
coldata[2,"id"] <- "SA10N"
coldata[18,"id"] <- "SA16D"
coldata[85,"id"] <- "SJ20D"
coldata[99,"id"] <- "SJ5O"

#remove duplicate samples in counts (keeping sample with most counts)
countssid <- subset(countssid, select = -c(SA12O.r_run2,SA12O,SA6O,SA6O.r_run2,SA9O))

#fix name issues in counts
colnames(countssid)[colnames(countssid) == "SA12O_run2"] <- "SA12O"
colnames(countssid)[colnames(countssid) == "SA6O_run2"] <- "SA6O"
colnames(countssid)[colnames(countssid) == "SA9O_run2"] <- "SA9O"

# Remove hybrid samples and juvenile samples from coldata
coldata_puread <- subset(coldata, age %in% "adult") %>%
  subset(!(admix %in% c("S2.D1")))

#remove individuals in mcavcounts that are not in the coldata list
countssid_common <- countssid[,colnames(countssid) %in% coldata_puread$id]

#remove individuals in coldata that are not in counts_common list
coldata_common <- coldata_puread[coldata_puread$id %in% colnames(countssid_common),]

#57 pure lineage, adult ssid samples with 2brad and tagseq data

# Remove genes with cumulative mean read count <5 
min.count <- 5
cum.geneCnt <- apply(countssid_common,1,mean)
table(cum.geneCnt >= min.count)
countssid_common <- countssid_common[cum.geneCnt >= min.count,]
#6610 genes > 5 counts

#Run deseq: needed for outlier detection
dds <- DESeqDataSetFromMatrix(countssid_common,
                              colData = coldata_common, 
                              design = formula(~ habitat + admix))

# Variance stabilizing transformation in order to minimize correlation between the variance and mean
Vsd <- vst(dds)
# creating normalized dataframe
vsd <- assay(Vsd)
# renames the column names
colnames(vsd) <- colnames(countssid_common)

# running outlier detection
errors=ExpressionSet(vsd, AnnotatedDataFrame(as.data.frame(colData(Vsd))))
arrayQualityMetrics(errors,intgroup=c("habitat","admix"), force=T)

#outliers: SA17N, SA20N, SA3O, SA4O, SA5O
countssid_common <- subset(countssid_common, select=-c(SA17N, SA20N, SA3O, SA4O, SA5O))
coldata_common <- subset(coldata_common, !(id %in% c("SA17N", "SA20N", "SA3O", "SA4O", "SA5O")))

#also removing only S2 individual: SA14N
countssid_common <- subset(countssid_common, select=-c(SA14N))
coldata_common <- subset(coldata_common, !(id %in% c("SA14N")))

#51 individuals remaining:
#habitat: 16 nearshore, 15 offshore, 20 deep
#lineage: 43 S1, 6 D1, 2 D2
#habitat x lineage: 16 nearshore S1, 15 offshore S1, 12 deep S1, 6 deep D1, 2 deep D2

# make the samples the row name for coldata
rownames(coldata_common) <- coldata_common[,1]
coldata_common$id <- NULL

#rerun array quality metrics without the crazy outlier to see if any other outliers exist in the data
dds <- DESeqDataSetFromMatrix(countssid_common,
                              colData = coldata_common, 
                              design = formula(~ habitat + admix))
# Variance stabilizing transformation in order to minimize correlation between the variance and mean
Vsd <- vst(dds)
# creating normalized dataframe
vsd <- assay(Vsd)
# renames the column names
colnames(vsd) <- colnames(countssid_common)

#save counts, coldata files
save(countssid_common, coldata_common, file = "ssid_files/counts_adult.Rdata")
save(dds, coldata_common, vsd, design, file="ssid_files/vsd_adult.RData")


