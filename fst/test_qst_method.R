library(data.table)
library(tidyverse)
library(dplyr)
library(reshape2)
library("doParallel")

load("~/Library/CloudStorage/Dropbox/Research/FLKeys_trans2_age/mcav_files/vsd_age.mcav.RData")

#create an incidence matrix
# Incidence matrix
incMat_a <- model.matrix(~ coldata_common$admix - 1)
colnames(incMat_a) <- gsub(pattern = 'ids\\$admix', replacement = "",
                           x = colnames(incMat_a))
rownames(incMat_a) <- rownames(coldata_common)
incMat_h <- model.matrix(~coldata_common$habitat - 1)
colnames(incMat_h) <- gsub(pattern = 'ids\\$habitat', replacement = "",
                           x = colnames(incMat_h))
rownames(incMat_h) <- rownames(coldata_common)

#load bams
bams <- read.table("bams_noclones", header = F)[,1]
bams <- gsub("_s2.trim.bt2.bam","",gsub("-r","",bams))

#load kinship matrix
kmat <- as.matrix(fread("newres",data.table = F))
kmat2 <- as.data.frame(subset(kmat, select = c("a","b","rab")))
kmat3 <- acast(kmat2, a~b, value.var='rab', fill=NULL, drop = F)
kmat3 <- cbind("0"=NA, kmat3)
kmat3 <- rbind(kmat3, "97"=NA)

#assign row and column names
rownames(kmat3) <- bams
colnames(kmat3) <- bams

#make the bottom diagonal mirror the top diagonal
makeSymm <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}
kmat4 <- makeSymm(kmat3)

#make NA = 1
kmat4[is.na(kmat4)] <- 1

#rename to match og script
K_ldak <- kmat4

#confusion here down:
K_ok <- K_ldak[match(unique(rownames(coldata_common)), rownames(K_ldak)),
               match(unique(rownames(coldata_common)), colnames(K_ldak))]
# Kw <- K_ok
# for (i in unique(pop_ids)) {
#   for (j in unique(pop_ids)) {
#     if (j != i) {
#       Kw[subset(ids, Population == i)$Genotype,
#          subset(ids, Population == j)$Genotype] <- 0
#     }
#   }
# }
# n <- ncol(Kw)
# KwNorm <- (n - 1) / sum((diag(n) - matrix(1, n, n) / n) * Kw) * Kw

# matrice d apparentement populations
Kb <- matrix(NA,
             nrow = length(levels(rownames(coldata_common))),
             ncol = length(levels(rownames(coldata_common))))
rownames(Kb) <- levels(rownames(coldata_common))
colnames(Kb) <- levels(rownames(coldata_common))
for (i in levels(rownames(coldata_common))) {
  for (j in levels(rownames(coldata_common))) {
    Kb[i, j] <- mean(K_ok[subset(coldata_common, habitat == i)$admix,
                          subset(coldata_common, habitat == j)$admix])
  }
}

n <- ncol(Kb)
KbNorm <- (n - 1) / sum((diag(n) - matrix(1, n, n) / n) * Kb) * Kb



## ----genotype model calculation v2, message=FALSE-------------------------------------------------------------------------------
#load counts data
load("~/Dropbox/Research/FLKeys_trans2_age/mcav_files/vsd_age.mcav.RData")

#assign a bunch of rqndom shit?
a <- proc.time()
numberOfThreads <- 24
cl <- makeCluster(numberOfThreads)
registerDoParallel(cl)
avancement <- 0

#extra random shit
genParFull <- c('(G_2_2_1_1+G_3_3_1_1)/(G_2_2_1_1+G_3_3_1_1+R_1_1)',
                'G_2_2_1_1/(G_2_2_1_1+2*G_3_3_1_1)')

#rename vsd to tcounts.trans (I think these are equivalent files counts that have been transformed)
tcounts.trans <- t(vsd)
genesAPasser <- 1:ncol(tcounts.trans)

#make vectors of coldata (what I think that those things are in colonnesASelectionner)
ids <- rownames(coldata_common)
age <- coldata_common$age
admix <- coldata_common$admix
habitat <- coldata_common$habiatat

colonnesASelectionner <- c("ids","age","admix","habitat")

#make rawCovsOk
rawCovsOk <- coldata_common
rawCovsOk$ids <- rownames(coldata_common)

if (all(rawCovsOk$Genotype_Bloc == rownames(tcounts.trans) %>% strsplit("_") %>%
        lapply(function(x) paste(x[2:3], collapse = "_")))) {
  data_breedR <- cbind(as.data.frame(rawCovsOk)[,colonnesASelectionner],
                       tcounts.trans)
} else {
  break()
}

data_breedR$Heure_prelev <- as.factor(format(round.POSIXt(as.POSIXct(
  data_breedR$Heure_prelev, format = "%H:%M", tz = "UTC"), units = "hours"),
  format = "%H:%M"))
data_breedR$Bloc <- as.factor(data_breedR$Bloc)
data_breedR <- data_breedR[order(data_breedR$Genotype_Bloc),]

varList <- foreach(j = genesAPasser,
                   .packages = c("breedR")) %dopar%
  {
    # avancement <- avancement + 1
    dataset <- data_breedR[, c(1:4, j + 4)]
    colnames(dataset)[c(1,5)] <- c("ids", "expr")
    message(j)
    # 
    #   dataset$expr <- tcounts.trans[, j]
    mod1 <- NULL
    try(mod1 <- suppressMessages(
      remlf90(fixed = expr ~ DatePrelev + Heure_prelev + Bloc, 
              generic = list(pop_ids = list(incMat_p, KbNorm),
                             genot_ids = list(incMat_g, KwNorm)),
              data = dataset,
              method = "em")),
      silent = TRUE)
    if (is.null(mod1)) { break() }
    blupGenotPlusPop <- incMat_gg %*% as.matrix(
      mod1$ranef$pop_ids[[1]]$value) + mod1$ranef$genot_ids[[1]]$value
    blupGenotPlusPop <- blupGenotPlusPop[,1]
    blupGenot <- mod1$ranef$genot_ids[[1]]$value
    names(blupGenot) <- names(blupGenotPlusPop)
    phenAjusted <- mod1$effects$pop_ids$effects[[1]]$incidence.matrix %*%
      as.matrix(mod1$ranef$pop_ids[[1]]$value) + 
      mod1$effects$genot_ids$effects[[1]]$incidence.matrix %*%
      as.matrix(mod1$ranef$genot_ids[[1]]$value) +
      residuals(mod1) + mean(dataset$expr)
    nPA <- rownames(phenAjusted)
    phenAjusted <- as.vector(phenAjusted)
    names(phenAjusted) <- nPA
    moyAjParGenot <- aggregate(phenAjusted ~ gsub("_[13]$", "", genot_id), mean,
                               data = dataset)
    genot <- moyAjParGenot[,1]
    moyAjParGenot <- moyAjParGenot[,2]
    names(moyAjParGenot) <- genot
    list(Variances = mod1$var[, "Estimated variances"],
         `BLUP Genot + Pop` = blupGenotPlusPop,
         `BLUP Genot` = blupGenot,
         AIC_wGenot = mod1$fit$AIC,
         `Adjusted Means Per Genotype` = moyAjParGenot,
         `Adjusted Phenotype` = phenAjusted)
  }
stopCluster(cl)
proc.time() - a



## ----nommage de la liste, cache=FALSE-------------------------------------------------------------------------------------------
names(varList) <- colnames(tcounts.trans)


## ----enlever les blups nuls-----------------------------------------------------------------------------------------------------
toRemove <- names(which(sapply(varList, function(x)
  sum(abs(x$`BLUP Genot + Pop`))) == 0))
varList <- varList[!names(varList) %in% toRemove]


## ----calcul des h2 et qst-------------------------------------------------------------------------------------------------------
h2EtQst <- as.data.frame(t(sapply(seq(along = varList), function(x){
  result = unname(varList[[x]]$Variances)
  c(h2 = sum(result[1:2])/sum(result[1:3]),
    qst = result[1]/(result[1] + 2*result[2]))
})))
par(mfrow = c(1,2))
hist(h2EtQst$h2)
hist(h2EtQst$qst)


