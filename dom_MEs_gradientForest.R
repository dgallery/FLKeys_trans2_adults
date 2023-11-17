library(ggplot2)
#install.packages("vegan")
library(vegan)
#install.packages("gradientForest", repos="http://R-Forge.R-project.org")
library(gradientForest)
#library(packfor)
library(viridis)
#library(DEseq2)
library(cowplot)


# environmental parameters (un-remark the species you want)
#meta=read.table("mcav_files/mc_enviro.meta",header=T)
meta=read.table("ssid_files/ss_enviro.meta",header=T)

# reading WGCNA results; MEs is what we need
#ll=load("mcav_files/networkdata_signed.adult.mcav.RData")
ll=load("ssid_files/networkdata_signed.adult.ssid.RData")
head(MEs)
# making sure metadata match the ME table
goods=row.names(MEs)[which(row.names(MEs) %in% meta$id)]
MEs=MEs[goods,]
row.names(meta)=meta$id
# leaving only habitat and lineage affiliation in metadata
meta=meta[goods,c(3,4)]

# oading vsd - we don't need it for ME analysis actually
ll=load("mcav_files/vsd_adults.mcav.RData")
#ll=load("vsd_adult.ssid.RData")
vsdt=t(vsd)
vsdt=vsdt[goods,]


# ---------------- dummifying metadata (turning charcter factors into 0-1)

covars.g=meta 
covs=c();ci=1
for (ci in 1:ncol(covars.g)) {
  if(is.factor(covars.g[,ci]) | is.integer(covars.g[,ci]) | is.character(covars.g[,ci])) { 
    co=as.factor(covars.g[,ci]) 
    covs=cbind(covs,model.matrix(~0+co)) 
  }  else { 
    co=covars.g[,ci]
    covs=cbind(covs,co)
    colnames(covs)[ncol(covs)]=colnames(covars.g)[ci]
  }
}	 
row.names(covs)=row.names(covars.g)
colnames(covs)=sub("^co","",colnames(covs))
covs=data.frame(covs)
covs

#----- RDA forward selection - skip that

# dev.off()
# ord.all=rda(MEs~.,covs)
# vsdist=dist(vsdt,method="manhattan")/ncol(vsdt)
# ord.all=capscale(vsdist~.,covs)
# plot(ord.all,scaling=1)
# r2adj.poly=RsquareAdj(ord.all)$adj.r.squared
# r2adj.poly
# 
# # forward selection
# fs.p=forward.sel(MEs,covs,adjR2thresh = r2adj.poly)
# fs.p$variables=factor(fs.p$variables,levels=(fs.p$variables[order(fs.p$R2)]))
# ggplot(fs.p,aes(variables,R2))+ geom_bar(stat="identity")+coord_flip()+theme_bw()
# 
# 
# capp=rda(MEs~.,covs[,fs.p[,2]])
# plot(capp,scaling=1)
# anova.cca(capp)
# capp.fit=data.frame(scores(capp,choices=c(1,2),scaling=1)$sites)

#---------------- gradient forest

# function to gradient forest PCs (all of them) of some ordination:
makeGF=function(ordination,environment,ntrees=500) {
  Y=data.frame(scores(ordination,scaling=1,choices=c(1:length(ordination$CA$eig)))$sites)
  X=environment
  nSites = dim(Y)[1]
  nSpecs = dim(Y)[2]
  lev = floor(log2(nSites * 0.368/2))
  gf = gradientForest(cbind(X, Y),
                      predictor.vars = colnames(X), response.vars = colnames(Y),
                      ntree = ntrees, transform = NULL, trace=T,
                      maxLevel = lev, corr.threshold = 0.25)
  return(gf)
}

# simpler function to gradient forest the raw dataset (without ordinating it)
makeGFa=function(Y,environment,ntrees=500) {
  X=environment
  nSites = dim(Y)[1]
  nSpecs = dim(Y)[2]
  lev = floor(log2(nSites * 0.368/2))
  gf = gradientForest(cbind(X, Y),
                      predictor.vars = colnames(X), response.vars = colnames(Y),
                      ntree = ntrees, transform = NULL, trace=T,
                      maxLevel = lev, corr.threshold = 0.25)
  return(gf)
}

# # computing ordination based on vsd - average manhattan distance (i.e. average log-fold change across genes)
# vsdist=dist(vsdt,method="manhattan")/ncol(vsdt)
# ord.all=capscale(vsdist~1)
# # gradient forest based on ordination 
# gf=makeGF(ord.all,covs,ntrees=1500)

# gradient forest straight on MEs
gf=makeGFa(MEs,covs,ntrees=1500)

# plotting bar chart of importances
ii=data.frame(importance(gf))
names(ii)="importance"
# reordering by decreasing importance
ii$var=factor(row.names(ii),levels=rev(row.names(ii)))

ggplot(ii[1:min(10,nrow(ii)),],aes(var,importance))+
  geom_bar(stat="identity")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_flip()+
#  scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
  xlab("")

