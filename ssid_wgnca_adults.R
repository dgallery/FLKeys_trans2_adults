#### wgcna for adult samples
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

load("ssid_files/counts_adult.Rdata")
load("ssid_files/vsd_adult.RData")

#--------- Run WGCNA --------------
#make datTraits and datExpr files
datTraits = coldata_common
datExpr = as.data.frame(t(vsd))

#add admix.habitat to datTraits
datTraits$admix.habitat <- paste0(datTraits$admix,datTraits$habitat, sep = "")

#make datTraits numeric
datTraits$D = ifelse(datTraits$habitat== "deep",1,0)
datTraits$N = ifelse(datTraits$habitat== "nearshore",1,0)
datTraits$O = ifelse(datTraits$habitat== "offshore",1,0)

datTraits$C1 = ifelse(datTraits$admix== "S1",1,0)
datTraits$C3 = ifelse(datTraits$admix== "D1",1,0)
datTraits$C4 = ifelse(datTraits$admix== "D2",1,0)

datTraits$S1N = ifelse(datTraits$admix.habitat=="S1nearshore",1,0)
datTraits$S1O = ifelse(datTraits$admix.habitat=="S1offshore",1,0)
datTraits$S1D = ifelse(datTraits$admix.habitat=="S1deep",1,0)
datTraits$D1D = ifelse(datTraits$admix.habitat=="D1deep",1,0)
datTraits$D2D = ifelse(datTraits$admix.habitat=="D2deep",1,0)

datTraits = datTraits %>%
  select(D,N,O,C1,C3,C4,S1N,S1O,S1D,D1D,D2D)
save(datTraits, datExpr, file = "ssid_files/wgcna_adults_input.Rdata")

# rename data frames to match pipeline
vsd.wg = vsd
datt=t(vsd)
traits=coldata_common
powers = c(seq(from = 2, to=26, by=2))
save(datt,traits,file="ssid_files/wgcnaData.adult.RData")

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType="signed")

par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

s.th=10
adjacency = adjacency(datExpr, power = s.th,type="signed");
TOM = TOMsimilarity(adjacency,TOMType="signed");
dissTOM = 1-TOM
save(TOM,dissTOM,file="ssid_files/TOM.adult.rdata")

# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30; 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average");

#merge modules
MEDissThres = 0.6 # in the first pass, set this to 0 - no merging (we want to see the module-traits heatmap first, then decide which modules are telling us the same story and better be merged)
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")  # on 2nd pass: does this cut height meet your merging goals? If not, reset MEDissThres and replot

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plotting the fabulous ridiculogram
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = FALSE, guideHang = 0.05,lwd=0.3)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# how many genes in each module?
table(moduleColors)

# Save module colors and labels for use in subsequent parts
save(MEs, geneTree, moduleLabels, moduleColors, file = "ssid_files/networkdata_signed.adult.RData")

#####Make the primary heatmap of module-trait relationships with the proper names#####
load("ssid_files/networkdata_signed.adult.ssid.RData")
load("ssid_files/wgcnaData.adult.RData")
load("ssid_files/wgcna_adults_input.Rdata")

#select the traits to graph and make the labels
#first heatmap = admix and site
Traits=datTraits[,c("C1", "C3", "C4","N", "O", "D")]
traitlabels <- c("Shallow 1\nEcomorph", "Deep 1\nEcomorph", "Deep 2\nEcomorph","Nearshore\nHabitat", "Offshore\nHabitat", "Deep\nHabitat")

nSamples = nrow(datt)
moduleTraitCor = cor(MEs, Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
head(MEs)

textMatrix = paste0(signif(moduleTraitCor, 2), "\n",
                    ifelse(moduleTraitPvalue > 0.05, "NS", 
                           ifelse(moduleTraitPvalue > 0.01, "*",
                                  ifelse(moduleTraitPvalue > 0.001, "**", "***"))))
dim(textMatrix) = dim(moduleTraitCor)

textMatrix2 = ifelse(moduleTraitPvalue > 0.05, "", 
                     ifelse(moduleTraitPvalue > 0.01, "*",
                            ifelse(moduleTraitPvalue > 0.001, "**", "***")))

heatlabels <- paste0(Hmisc::capitalize(sub("ME","",row.names(moduleTraitCor))), '\n(n=', as.vector(table(moduleColors)[sub("ME","",row.names(moduleTraitCor))]), ')')
par(mar=c(6.5,7.5,2,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = traitlabels,
               xLabelsAngle = 60,
               #xLabelsAdj = c(0.5,0),
               #    yLabels = names(MEs),
               yLabels = heatlabels,
               ySymbols = row.names(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix2,
               setStdMargins = FALSE,
               cex.text = 1.3,
               cex.lab.x = 1.3,
               cex.lab.y = 1.3,
               textAdj = c(0.5, 0.6),
               cex.legendLabel = 1.2,
               zlim = c(-1,1))

#plot only single module
# keeps="MEbrown"
# df=as.matrix(moduleTraitCor[row.names(moduleTraitCor) %in% keeps,5:7])
# df=t(df)
# m.sizes=c("Brown\n(n = 3987)")
# 
# keeps="6"
# tmatrix=as.matrix(textMatrix[row.names(textMatrix) %in% keeps,5:7])
# tmatrix=t(tmatrix)
# trait_lab=traitlabels[5:7]
# 
# labeledHeatmap(Matrix = df,
#                xLabels = trait_lab,
#                xLabelsAngle = 60,
#                yLabels = m.sizes,
#                ySymbols = row.names(df),
#                colorLabels = FALSE,
#                colors = blueWhiteRed(50),
#                textMatrix = tmatrix,
#                setStdMargins = FALSE,
#                cex.text = 1.3,
#                cex.lab.x = 1,
#                cex.lab.y = 1,
#                textAdj = c(0.5, 0.6),
#                cex.legendLabel = 0.6,
#                zlim = c(-1,1))
# 

# Plot Heatmap of admix/habitat combos
Traits2=datTraits[,c("S1N", "S1O", "S1D", "D1D", "D2D")]
traitlabels <- c("Shallow 1 Ecomorph\nNearshore Habitat", "Shallow 1 Ecomorph\nOffshore Habitat", "Shallow 1 Ecomorph\nDeep Habitat",
                 "Deep 1 Ecomorph\nDeep Habitat", "Deep 2 Ecomorph\nDeep Habitat")

nSamples = nrow(datt)
moduleTraitCor = cor(MEs, Traits2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
head(MEs)

textMatrix = paste0(signif(moduleTraitCor, 2), "\n",
                    ifelse(moduleTraitPvalue > 0.05, "NS", 
                           ifelse(moduleTraitPvalue > 0.01, "*",
                                  ifelse(moduleTraitPvalue > 0.001, "**", "***"))))
dim(textMatrix) = dim(moduleTraitCor)
textMatrix2 = ifelse(moduleTraitPvalue > 0.05, "", 
                     ifelse(moduleTraitPvalue > 0.01, "*",
                            ifelse(moduleTraitPvalue > 0.001, "**", "***")))


heatlabels <- paste0(Hmisc::capitalize(sub("ME","",row.names(moduleTraitCor))), '\n(n=', as.vector(table(moduleColors)[sub("ME","",row.names(moduleTraitCor))]), ')')
par(mar=c(9,7.5,2,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = traitlabels,
               xLabelsAngle = 60,
               yLabels = heatlabels,
               ySymbols = row.names(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix2,
               setStdMargins = FALSE,
               cex.text = 1,
               cex.lab.x = 1,
               cex.lab.y = 1,
               textAdj = c(0.5, 0.6),
               cex.legendLabel = 1.2,
               zlim = c(-1,1))

#-------barplots of modules ordered by site------
load("ssid_files/networkdata_signed.adult.RData")
load("ssid_files/wgcnaData.adult.RData")
load("ssid_files/wgcna_adults_input.Rdata")

datTraits2 <- datTraits[order(datTraits$D,datTraits$O,datTraits$N),]

which.module="grey60"
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
datME <- datME[match(row.names(datTraits2), rownames(datME)),]

ME=datME[, paste("ME",which.module, sep="")]
par(mar=c(4, 4.2, 4, 1))
barplot(ME, col=which.module, main="Grey60 Module", cex.main=1,
        ylab="Eigengene expression",xlab="",
        mgp=c(3,1,0), names.arg = rownames(datME), las =2, cex.names = .7)
mtext("Samples ordered by habitat", side = 1, line =3)

#barplots ordered by admix
datTraits3 <- datTraits[order(datTraits$C4,datTraits$C3,datTraits$C1),]

which.module="grey60"
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
datME <- datME[match(row.names(datTraits3), rownames(datME)),]

ME=datME[, paste("ME",which.module, sep="")]
par(mar=c(4, 4.2, 4, 1))
barplot(ME, col=which.module, main="Grey60 module", cex.main=1,
        ylab="Eigengene expression",xlab="",
        mgp=c(3,1,0), names.arg = rownames(datME), las =2, cex.names = 0.7)
mtext("Samples ordered by ecomorph", side = 1, line =3)

######### eigengene boxplots ##########
# y-axis = eigengene expression
# x-axis = habitat/lineage
load("ssid_files/counts_adult.Rdata")
load("ssid_files/networkdata_signed.adult.RData")

grey60_df <- merge(MEs, coldata_common, by = "row.names")
grey60_df <- na.omit(grey60_df)
##order by lineage
grey60_df$admix <- factor(grey60_df$admix, levels = c("S1","D1","D2"))
grey60_df <- na.omit(grey60_df)
ggplot(grey60_df, aes(x = admix, y = MEgrey60, color = admix))+
  geom_boxplot(width = 0.5, outlier.shape = NA)+
  theme_bw()+
  theme(panel.grid = element_blank(),  axis.text.x = element_blank())+
  ggforce::geom_sina(maxwidth = 0.3)+
  scale_color_manual(name = "Ecomorph", labels = c("Shallow 1","Deep 1", "Deep 2"), values = c('palegreen1','skyblue','dodgerblue4')) +
  scale_y_continuous(name=expression("Grey60 Module Eigengene"))+
  scale_x_discrete(name=expression("Ecomorph"))

#order by habitat
grey60_df$habitat <- factor(grey60_df$habitat, levels = c("nearshore","offshore","deep"))
ggplot(grey60_df, aes(x = habitat, y = MEgrey60, color = habitat))+
  geom_boxplot(width = 0.5, outlier.shape = NA)+
  theme_bw()+
  theme(panel.grid = element_blank(),  axis.text.x = element_blank())+
  ggforce::geom_sina(maxwidth = 0.3)+
  scale_color_manual(name = "Habitat", labels = c("Nearshore","Offshore", "Deep"), values = c('#004b23','#8b0a50','#002a52')) +
  scale_y_continuous(name=expression("Grey60 Module Eigengene"))+
  scale_x_discrete(name=expression("Habitat"))

#################
# saving  modules for GOMWU analysis (two-parts: Fisher test, MWU test within-module)
load("ssid_files/wgcnaData.adult.RData")
load("ssid_files/networkdata_signed.adult.RData")

vsd=t(datt)

# calculating modul memberships for all genes for all modules
allkME =as.data.frame(signedKME(datt, MEs)) 
names(allkME)=gsub("kME","",names(allkME))

whichModule="grey60" #make one for each module
table(moduleColors==whichModule) # how many genes are in the modules
vsd.wg=t(datt) 

# Saving data for Fisher-MWU combo test (GO_MWU - full version for between species comparison)
inModuleBinary=as.numeric(moduleColors==whichModule)
combo=data.frame("gene"=row.names(vsd.wg),"Fish_kME"=allkME[,whichModule])
write.csv(combo,file=paste("ssid_go/",whichModule,"_full.csv",sep=""),row.names=F,quote=F)

# Saving data for Fisher-MWU combo test (GO_MWU - normal version for within species comparison)
combo=data.frame("gene"=row.names(vsd.wg),"Fish_kME"=allkME[,whichModule]*inModuleBinary)
write.csv(combo,file=paste("ssid_go/",whichModule,".csv",sep=""),row.names=F,quote=F)
