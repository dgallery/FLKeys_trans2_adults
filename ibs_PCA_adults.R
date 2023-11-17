######## IBS #########
#~~~~~~~~~~~~~~~~~~~~~~~ 2bRAD ~~~~~~~~~~~~~~~~~~~~~~~~#

# 2bRAD_bothsets ----------------------------------------------------------

#both datasets to ID clones, wrong species, and assign clusters
library(vegan)
library(ggplot2)
library(dplyr)

spp <- 'ssid' #choose species
npops <- 4 #choose number of populations

dir <- paste0('popgen/',spp,"_both_transects/") # change this to where your scp'd files are
bams <- read.table(paste0(dir, "bams"))[,1] # list of bam files
if(spp=='ssid'){
  prefix <- 'ss'
  bams <- gsub(".nosymbio.fastq.bam","",bams)
} else if (spp=='mcav'){
  prefix <- 'mc'
  bams <- gsub(".bt2.bam","", bams)
}

ibsMat <- as.matrix(read.table(paste0(dir,prefix,"02ball.ibsMat")))
dimnames(ibsMat) <- list(bams,bams)

hc <- hclust(as.dist(ibsMat),"ave")
plot(hc,cex=0.7) # clustering of samples by IBS (great to detect clones or closely related individuals)

#mcav clones: 
if(spp=='mcav'){
  clones <- c("MJ6D_s2.trim", "MA6N_s2.trim", "MCDA20c-1.fastq", "MCDA20c-3.fastq", "MCOA11.fastq", 
              "MCOA14.fastq", "MCNJ5.fastq", "MCOA9.fastq","MCOA10c.fastq", "MCDJ19c-2.fastq", 
              "MCDJ19c-3.fastq", "MCNJ15c.fastq", "MCOJ5c.fastq", "MCOA15.fastq", 
              "MCOJ15.fastq", "MCNA16.fastq", "MCNA17.fastq", "MCNA10.fastq", "MCOA15.fastq")
} else if(spp == 'ssid'){
  clones <- c("SSNA11","SSNA1","SSDJ19c-2", "SSDJ19c-3","SSOJ20c-1", "SSOJ20c-2","SJ2Ob_s2",
              "SSNA8","SA10N_s2","SA16D_s2","SJ17N-r_s2","SA11N_s2","SJ8Ob_s2","SSNA10","SSNA9",
              "SSNJ20c-1", "SSNJ20c-2", "SSNJ20c-3","SJ1N_s2","SJ18N_s2","SJ19N_s2","SJ12N_s2","SJ13N_s2","SSNJ10",
              "SSDJ20c-1","SSDA20c-3", "SSDA20c-1")
}
#cross transect clones: removing the set 1 individuals, but leaving the set 2 individuals
#still have no idea HOW we have cross transect clones so that's fun

#remove clones from bams list:
goods <- which(!(bams %in% clones))
clonn <- which((bams %in% clones))
ibsMat <- ibsMat[goods,goods]
bams <- bams[goods]
if(spp=='ssid'){
  bamlist <- data.frame(paste(bams,".nosymbio.fastq.bam",sep=""))
} else if (spp=='mcav'){
  bamlist <- data.frame(paste(bams,".bt2.bam",sep=""))
}
write.table(bamlist,file=paste0(dir,"bams_noclones"),row.names=F,col.names=F,quote=F) # use this bam list to rerun ANGSD for ADMIXTURE

#---------------
#after rerunning ANGSD for ADMIXTURE
#clear environment before redoing with new IBS matrix without clones
spp <- 'mcav' #choose species
npops <- 4 #choose number of populations

dir <- paste0('popgen/',spp,"_both_transects/") # change this to where your scp'd files are
bams <- read.table(paste0(dir, "bams_noclones"))[,1] # list of bam files
if(spp=='ssid'){
  prefix <- 'ss'
  bams <- gsub(".nosymbio.fastq.bam","",bams)
} else if (spp=='mcav'){
  prefix <- 'mc'
  bams <- gsub(".trim.bt2.bam","", gsub(".fastq.bt2.bam","",bams))
}

ibsMat <- as.matrix(read.table(paste0(dir,spp,"2ball.ibsMat")))
dimnames(ibsMat) <- list(bams,bams)

hc <- hclust(as.dist(ibsMat),"ave")
plot(hc,cex=0.7) 

save(bams,ibsMat,file = paste0(dir,spp,"_2ball_ibs_bams.rdata"))

#--------get *clusters.rdata file with admixture.R to make metadata files and plot PCAs below------

# Make metadata and environment files ----------------------------------------------------------
spp <- 'mcav' #choose species
npops <- 4 #choose number of populations
if(spp=='ssid'){
  prefix <- 'ss'
} else if (spp=='mcav'){
  prefix <- 'mc'
}

load(paste0("popgen/",spp,"_both_transects/",spp,"_2ball_ibs_bams.rdata"))

adults=grep("A",bams)
juveniles=grep("J",bams)
near=grep("N",bams)
off=grep("O",bams)
deep=grep("D",bams)

# factors:
site=rep("deep",length(bams))
site[off]="offshore"
site[near]="nearshore"
aj=rep("adult",length(bams))
aj[juveniles]="juvenile"

load(paste0('popgen/',spp,"_both_transects/",prefix,"_k4.qopt_clusters.RData"))

set=rep(1,length(bams))
set2=grep("_s2",bams)
set[set2]=2
set=as.factor(set)

meta <- data.frame(id = bams, age = aj, habitat = site, admix = cluster.admix, transect = set)
if(spp=='mcav'){
  meta$id <- gsub("_s2","",meta$id)
} else if(spp=='ssid'){
  meta$id <- gsub("_s2","",meta$id)
}

#rename admix to match JP's paper:
if(spp=='mcav'){
  mc_meta <- meta %>%
  mutate(admix = recode(admix, "1"="D3", "2"="O2", "3"="N1", "4"="D4", "1.2"="OD2.3", 
                        "2.3"="NO1.2", "2.4"="OD2.4"))
}else if(spp=='ssid'){
  ss_meta <- meta %>%
    mutate(admix = recode(admix, "1"="D1","2"="D2","3"="S1","4"="S2","1.4"="S2.D1","1.3"="S1.D1","3.4"="S1.S2"))###fix these
}

#save file
if(spp=='mcav'){
  write.table(mc_meta,file="mcav_files/mcav.metadata.2ball",row.names=F,col.names=F,quote=F)
  save(mc_meta,file = "popgen/mcav_both_transects/mc_metadata_2ball.rdata")
}else if(spp=='ssid'){
  write.table(ss_meta,file="ssid_files/ssid.metadata.2ball",row.names=F,col.names=F,quote=F)
  save(ss_meta,file = "popgen/ssid_both_transects/ss_metadata_2ball.rdata")
}

#-----overlapping PCA plots of trans 1 and trans 2------
library(vegan)
library(ggplot2)

spp <- 'mcav' #choose species
npops <- 4 #choose number of populations
if(spp=='ssid'){
  prefix <- 'ss'
} else if (spp=='mcav'){
  prefix <- 'mc'
}

#load files:
if(spp == 'mcav'){
  load("popgen/mcav_both_transects/mcav_2ball_ibs_bams.rdata")
  load("popgen/mcav_both_transects/mc_metadata_2ball.rdata")
} else if (spp == 'ssid'){
  load("popgen/ssid_both_transects/ssid_2ball_ibs_bams.rdata")
  load("popgen/ssid_both_transects/ss_metadata_2ball.rdata")
}

#remove juveniles from bams and ibsmat
bams_ad <- bams[grep("A",bams)]
ibsMat_ad <- ibsMat[bams_ad,bams_ad]
if(spp=='mcav'){
  mc_meta_ad <- mc_meta[- grep("juvenile",mc_meta$age),]
}else if(spp=='ssid'){
  ss_meta_ad <- ss_meta[- grep("juvenile",ss_meta$age),]
}

#run capscale
pp0_ad=capscale(ibsMat_ad~1)

#make transect factors
set=rep(1,length(bams_ad))
set2=grep("_s2",bams_ad)
set[set2]=2
set=as.factor(set)

#set which transect to run
transect2plot=2
if(transect2plot==1){
  set <- factor(set, levels=c('1','2'))
} else if(transect2plot==2){
  set <- factor(set, levels=c('2','1'))
}

#set which mds to plot
axes2plot=c(1,2)
axes=1.2

if(spp=='mcav'){
  pp0.admix <- data.frame(id=mc_meta_ad$id, age=mc_meta_ad$age, habitat=mc_meta_ad$habitat, transect=set, admix=mc_meta_ad$admix, pp0_ad$CA$u)
  admix.spider <- merge(pp0.admix, aggregate(cbind(mean.x=pp0_ad$CA$u[,axes2plot[1]], mean.y=pp0_ad$CA$u[,axes2plot[2]])~admix, pp0.admix, mean), by='admix')
  eigvals <- summary(eigenvals(pp0_ad))
}else if(spp=='ssid'){
  pp0.admix <- data.frame(id=ss_meta_ad$id, age=ss_meta_ad$age, habitat=ss_meta_ad$habitat, transect=set, admix=ss_meta_ad$admix, pp0_ad$CA$u)
  admix.spider <- merge(pp0.admix, aggregate(cbind(mean.x=pp0_ad$CA$u[,axes2plot[1]], mean.y=pp0_ad$CA$u[,axes2plot[2]])~admix, pp0.admix, mean), by='admix')
  eigvals <- summary(eigenvals(pp0_ad))
}
#reorder admix.spider by habitat
if(spp == 'mcav'){
  admix.spider$admix <- factor(admix.spider$admix, levels = c("N1","O2","D3","D4","NO1.2","OD2.3","OD2.4"))
  admix.spider$transect <- factor(admix.spider$transect, levels = c(2,1))
} else if (spp == 'ssid'){
  admix.spider$admix <- factor(admix.spider$admix, levels = c("S1","S2","D1","D2")) #no hybrids in adult ssids
  admix.spider$transect <- factor(admix.spider$transect, levels = c(2,1))
}

#ggplot
# if(spp=='mcav' && transect2plot=="1" && axes=="1.2"){
#   mc_trans1_pc12 <- ggplot(data=admix.spider)+
#     theme_bw()+
#     theme(panel.grid = element_blank())+
#     labs(x = paste0('MDS', axes2plot[1], ' (', round(eigvals[2,axes2plot[1]]*100,1), '%)'), y = paste0('MDS', axes2plot[2], ' (', round(eigvals[2,axes2plot[2]]*100,1), '%)'))+
#     geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], yend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], alpha='transect'), col='grey70')+
#     geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape='age', fill='admix', col='transect', alpha='transect'), size = 5)+
#     scale_color_manual(values = c('black','white'))+
#     scale_alpha_manual(values = c(1,0.2)) +
#     scale_fill_manual(values = c('palegreen1', 'plum1','skyblue1','dodgerblue4','black','grey40','grey80'))+
#     scale_shape_manual(values = c(21,24))+
#     #coord_fixed(xlim = c(-0.15,0.3))+ #ssid
#     #coord_fixed(xlim = c(-0.22,0.1))+  #mcav pc1,2
#     theme(legend.position = "none", axis.title = element_text(size =16),axis.text = element_text(size = 14))
#     #guides(fill=guide_legend(override.aes = list(colour=c('palegreen1', 'plum1','skyblue1','dodgerblue4'))), shape = "none", alpha = "none", col = "none")
#   mc_trans1_pc12
# } else if(spp=='mcav' && transect2plot=="2" && axes=="1.2"){
#   mc_trans2_pc12 <- ggplot(data=admix.spider)+
#     theme_bw()+
#     theme(panel.grid = element_blank())+
#     labs(x = paste0('MDS', axes2plot[1], ' (', round(eigvals[2,axes2plot[1]]*100,1), '%)'), y = paste0('MDS', axes2plot[2], ' (', round(eigvals[2,axes2plot[2]]*100,1), '%)'))+
#     geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], yend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], alpha='transect'), col='grey70')+
#     geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape='age', fill='admix', col='transect', alpha='transect'), size = 5)+
#     scale_color_manual(values = c('black','white'))+
#     scale_alpha_manual(values = c(1,0.2)) +
#     scale_fill_manual(values = c('palegreen1', 'plum1','skyblue1','dodgerblue4','black','grey40','grey80'))+
#     scale_shape_manual(values = c(21,24))+
#     #coord_fixed(xlim = c(-0.15,0.3))+ #ssid
#     #coord_fixed(xlim = c(-0.22,0.1))+  #mcav pc1,2
#     theme(legend.position = "none", axis.title = element_text(size =16),axis.text = element_text(size = 14))
#   #guides(fill=guide_legend(override.aes = list(colour=c('palegreen1', 'plum1','skyblue1','dodgerblue4'))), shape = "none", alpha = "none", col = "none")
#   mc_trans2_pc12
# }else if(spp=='mcav' && transect2plot=="1" && axes=="3.2"){
#   mc_trans1_pc32 <- ggplot(data=admix.spider)+
#     theme_bw()+
#     theme(panel.grid = element_blank())+
#     labs(x = paste0('MDS', axes2plot[1], ' (', round(eigvals[2,axes2plot[1]]*100,1), '%)'), y = paste0('MDS', axes2plot[2], ' (', round(eigvals[2,axes2plot[2]]*100,1), '%)'))+
#     geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], yend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], alpha='transect'), col='grey70')+
#     geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape='age', fill='admix', col='transect', alpha='transect'), size = 5)+
#     scale_color_manual(values = c('black','white'))+
#     scale_alpha_manual(values = c(1,0.2)) +
#     scale_fill_manual(values = c('palegreen1', 'plum1','skyblue1','dodgerblue4','black','grey40','grey80'))+
#     scale_shape_manual(values = c(21,24))+
#     #coord_fixed(xlim = c(-0.15,0.3))+ #ssid
#     #coord_fixed(xlim = c(-0.22,0.1))+  #mcav pc1,2
#     theme(legend.position = "none", axis.title = element_text(size =16),axis.text = element_text(size = 14))
#   #guides(fill=guide_legend(override.aes = list(colour=c('palegreen1', 'plum1','skyblue1','dodgerblue4'))), shape = "none", alpha = "none", col = "none")
#   mc_trans1_pc32
# }else if(spp=='mcav' && transect2plot=="2" && axes=="3.2"){
#   mc_trans2_pc32 <- ggplot(data=admix.spider)+
#     theme_bw()+
#     theme(panel.grid = element_blank())+
#     labs(x = paste0('MDS', axes2plot[1], ' (', round(eigvals[2,axes2plot[1]]*100,1), '%)'), y = paste0('MDS', axes2plot[2], ' (', round(eigvals[2,axes2plot[2]]*100,1), '%)'))+
#     geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], yend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], alpha='transect'), col='grey70')+
#     geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape='age', fill='admix', col='transect', alpha='transect'), size = 5)+
#     scale_color_manual(values = c('black','white'))+
#     scale_alpha_manual(values = c(1,0.2)) +
#     scale_fill_manual(values = c('palegreen1', 'plum1','skyblue1','dodgerblue4','black','grey40','grey80'))+
#     scale_shape_manual(values = c(21,24))+
#     #coord_fixed(xlim = c(-0.15,0.3))+ #ssid
#     #coord_fixed(xlim = c(-0.22,0.1))+  #mcav pc1,2
#     theme(legend.position = "none", axis.title = element_text(size =16),axis.text = element_text(size = 14))
#   #guides(fill=guide_legend(override.aes = list(colour=c('palegreen1', 'plum1','skyblue1','dodgerblue4'))), shape = "none", alpha = "none", col = "none")
#   mc_trans2_pc32
# }else if(spp=='ssid' && transect2plot=="1" && axes=="1.2"){
#   ss_trans1_pc12 <- ggplot(data=admix.spider)+
#     theme_bw()+
#     theme(panel.grid = element_blank())+
#     labs(x = paste0('MDS', axes2plot[1], ' (', round(eigvals[2,axes2plot[1]]*100,1), '%)'), y = paste0('MDS', axes2plot[2], ' (', round(eigvals[2,axes2plot[2]]*100,1), '%)'))+
#     geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], yend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], alpha='transect'), col='grey70')+
#     geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape='age', fill='admix', col='transect', alpha='transect'), size = 5)+
#     scale_color_manual(values = c('black','white'))+
#     scale_alpha_manual(values = c(1,0.2)) +
#     scale_fill_manual(values = c('palegreen1', 'plum1','skyblue1','dodgerblue4'))+
#     scale_shape_manual(values = c(21,24))+
#     #coord_fixed(xlim = c(-0.15,0.3))+ #ssid
#     #coord_fixed(xlim = c(-0.22,0.1))+  #ssid pc1,2
#     theme(legend.position = "none", axis.title = element_text(size =16),axis.text = element_text(size = 14))
#   #guides(fill=guide_legend(override.aes = list(colour=c('palegreen1', 'plum1','skyblue1','dodgerblue4'))), shape = "none", alpha = "none", col = "none")
#   ss_trans1_pc12
# } else if(spp=='ssid' && transect2plot=="2" && axes=="1.2"){
#   ss_trans2_pc12 <- ggplot(data=admix.spider)+
#     theme_bw()+
#     theme(panel.grid = element_blank())+
#     labs(x = paste0('MDS', axes2plot[1], ' (', round(eigvals[2,axes2plot[1]]*100,1), '%)'), y = paste0('MDS', axes2plot[2], ' (', round(eigvals[2,axes2plot[2]]*100,1), '%)'))+
#     geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], yend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], alpha='transect'), col='grey70')+
#     geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape='age', fill='admix', col='transect', alpha='transect'), size = 5)+
#     scale_color_manual(values = c('black','white'))+
#     scale_alpha_manual(values = c(1,0.2)) +
#     scale_fill_manual(values = c('palegreen1', 'plum1','skyblue1','dodgerblue4'))+
#     scale_shape_manual(values = c(21,24))+
#     #coord_fixed(xlim = c(-0.15,0.3))+ #ssid
#     #coord_fixed(xlim = c(-0.22,0.1))+  #ssid pc1,2
#     theme(legend.position = "none", axis.title = element_text(size =16),axis.text = element_text(size = 14))
#   #guides(fill=guide_legend(override.aes = list(colour=c('palegreen1', 'plum1','skyblue1','dodgerblue4'))), shape = "none", alpha = "none", col = "none")
#   ss_trans2_pc12
# }else if(spp=='ssid' && transect2plot=="1" && axes=="3.2"){
#   ss_trans1_pc32 <- ggplot(data=admix.spider)+
#     theme_bw()+
#     theme(panel.grid = element_blank())+
#     labs(x = paste0('MDS', axes2plot[1], ' (', round(eigvals[2,axes2plot[1]]*100,1), '%)'), y = paste0('MDS', axes2plot[2], ' (', round(eigvals[2,axes2plot[2]]*100,1), '%)'))+
#     geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], yend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], alpha='transect'), col='grey70')+
#     geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape='age', fill='admix', col='transect', alpha='transect'), size = 5)+
#     scale_color_manual(values = c('black','white'))+
#     scale_alpha_manual(values = c(1,0.2)) +
#     scale_fill_manual(values = c('palegreen1', 'plum1','skyblue1','dodgerblue4'))+
#     scale_shape_manual(values = c(21,24))+
#     #coord_fixed(xlim = c(-0.15,0.3))+ #ssid
#     #coord_fixed(xlim = c(-0.22,0.1))+  #ssid pc1,2
#     theme(legend.position = "none", axis.title = element_text(size =16),axis.text = element_text(size = 14))
#   #guides(fill=guide_legend(override.aes = list(colour=c('palegreen1', 'plum1','skyblue1','dodgerblue4'))), shape = "none", alpha = "none", col = "none")
#   ss_trans1_pc32
# }else if(spp=='ssid' && transect2plot=="2" && axes=="3.2"){
#   ss_trans2_pc32 <- ggplot(data=admix.spider)+
#     theme_bw()+
#     theme(panel.grid = element_blank())+
#     labs(x = paste0('MDS', axes2plot[1], ' (', round(eigvals[2,axes2plot[1]]*100,1), '%)'), y = paste0('MDS', axes2plot[2], ' (', round(eigvals[2,axes2plot[2]]*100,1), '%)'))+
#     geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], yend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], alpha='transect'), col='grey70')+
#     geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape='age', fill='admix', col='transect', alpha='transect'), size = 5)+
#     scale_color_manual(values = c('black','white'))+
#     scale_alpha_manual(values = c(1,0.2)) +
#     scale_fill_manual(values = c('palegreen1', 'plum1','skyblue1','dodgerblue4'))+
#     scale_shape_manual(values = c(21,24))+
#     #coord_fixed(xlim = c(-0.15,0.3))+ #ssid
#     #coord_fixed(xlim = c(-0.22,0.1))+  #mcav pc1,2
#     theme(legend.position = "none", axis.title = element_text(size =16),axis.text = element_text(size = 14))
#   #guides(fill=guide_legend(override.aes = list(colour=c('palegreen1', 'plum1','skyblue1','dodgerblue4'))), shape = "none", alpha = "none", col = "none")
#   ss_trans2_pc32
# }
# 
# #legends
# if(spp=='mcav'){
#   legend_pca_mc <- ggplot(data=admix.spider)+
#     geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], fill='admix'), size = 5)+
#     lims(x= c(0,0), y= c(0,0))+
#     theme_void()+
#     theme(legend.position = c(0.4,0.5), legend.text = element_text(size =  12),
#           legend.title = element_text(size = 15))+
#     scale_fill_manual(values=c('palegreen1', 'plum1','skyblue1','dodgerblue4','black','grey40','grey80'),
#                       labels = c("Nearshore","Offshore","Deep1","Deep2","Nearshore-Offshore","Offshore-Deep1","Offshore-Deep2"))+
#     guides(fill=guide_legend(title = "Lineage",override.aes = list(shape = 21)))
#   legend_pca_mc
# }else if(spp=='ssid'){
#   legend_pca_ss <- ggplot(data=admix.spider)+
#     geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], fill='admix'), size = 5)+
#     lims(x= c(0,0), y= c(0,0))+
#     theme_void()+
#     theme(legend.position = c(0.2,0.5), legend.text = element_text(size =  12),
#           legend.title = element_text(size = 15))+
#     scale_fill_manual(values=c('palegreen1', 'plum1','skyblue1','dodgerblue4','black','grey40','grey80'),
#                       labels = c("Shallow1","Shallow2","Deep1","Deep2"))+
#     guides(fill=guide_legend(title = "Ecomorph", override.aes = list(shape = 21)))
#   legend_pca_ss
# }
# 
# if(spp == 'mcav'){
#   save(mc_trans1_pc12,mc_trans1_pc32,mc_trans2_pc12,mc_trans2_pc32,legend_pca_mc,legend_pca_mc, file = "popgen/mcav_bothsets_pca.rdata") 
# }else if(spp == 'ssid'){
#   save(ss_trans1_pc12,ss_trans1_pc32,ss_trans2_pc12,ss_trans2_pc32,legend_pca_ss, file = "popgen/ssid_bothsets_pca.rdata") 
# }

##### overlaping pca plots #####


mc_bothtrans_pc12 <- ggplot(data=admix.spider)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    labs(x = paste0('MDS', axes2plot[1], ' (', round(eigvals[2,axes2plot[1]]*100,1), '%)'), y = paste0('MDS', axes2plot[2], ' (', round(eigvals[2,axes2plot[2]]*100,1), '%)'))+
    geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], yend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], alpha='transect'), col='grey70')+
    geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape='transect', fill='admix', col='transect', alpha='transect'), size = 5)+
    scale_color_manual(values = c('black','white'))+
    scale_alpha_manual(values = c(1,0.7)) +
    scale_fill_manual(values = c('palegreen1', 'plum1','skyblue1','dodgerblue4','black','grey40','grey80'))+
    scale_shape_manual(values = c(21,24))+
    #coord_fixed(xlim = c(-0.15,0.3))+ #ssid
    #coord_fixed(xlim = c(-0.22,0.1))+  #mcav pc1,2
    theme(legend.position = "none", axis.title = element_text(size =16),axis.text = element_text(size = 14))
  #guides(fill=guide_legend(override.aes = list(colour=c('palegreen1', 'plum1','skyblue1','dodgerblue4'))), shape = "none", alpha = "none", col = "none")
mc_bothtrans_pc12
ss_bothtrans_pc12 <- ggplot(data=admix.spider)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    labs(x = paste0('MDS', axes2plot[1], ' (', round(eigvals[2,axes2plot[1]]*100,1), '%)'), y = paste0('MDS', axes2plot[2], ' (', round(eigvals[2,axes2plot[2]]*100,1), '%)'))+
    geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], yend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], alpha='transect'), col='grey70')+
    geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape='transect', fill='admix', col='transect', alpha='transect'), size = 5)+
    scale_color_manual(values = c('black','white'))+
    scale_alpha_manual(values = c(1,0.7)) +
    scale_fill_manual(values = c('palegreen1', 'plum1','skyblue1','dodgerblue4'))+
    scale_shape_manual(values = c(21,24))+
    #coord_fixed(xlim = c(-0.15,0.3))+ #ssid
    #coord_fixed(xlim = c(-0.22,0.1))+  #ssid pc1,2
    theme(legend.position = "none", axis.title = element_text(size =16),axis.text = element_text(size = 14))
  #guides(fill=guide_legend(override.aes = list(colour=c('palegreen1', 'plum1','skyblue1','dodgerblue4'))), shape = "none", alpha = "none", col = "none")
ss_bothtrans_pc12


#legends
legend_pca_mc_1 <- ggplot(data=admix.spider)+
    geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], fill='admix'), size = 5)+
    lims(x= c(0,0), y= c(0,0))+
    theme_void()+
    theme(legend.position = c(0.4,0.5), legend.text = element_text(size =  12),
          legend.title = element_text(size = 15))+
    scale_fill_manual(values=c('palegreen1', 'plum1','skyblue1','dodgerblue4','black','grey40','grey80'),
                      labels = c("Nearshore","Offshore","Deep1","Deep2","Nearshore-Offshore","Offshore-Deep1","Offshore-Deep2"))+
    guides(fill=guide_legend(title = "Lineage",override.aes = list(shape = 21)))
legend_pca_mc_1
legend_pca_mc_2 <- ggplot(data=admix.spider)+
    geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape = "transect"), size = 5)+
    lims(x= c(0,0), y= c(0,0))+
    theme_void()+
    theme(legend.position = c(0.4,0.5), legend.text = element_text(size =  12),
          legend.title = element_text(size = 15))+
    scale_shape_manual(values=c(16,17),labels = c("West","East"))+
    guides(shape=guide_legend(title = "Transect"))
legend_pca_mc_2
legend_pca_ss_1 <- ggplot(data=admix.spider)+
    geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], fill='admix'), size = 5)+
    lims(x= c(0,0), y= c(0,0))+
    theme_void()+
    theme(legend.position = c(0.2,0.5), legend.text = element_text(size =  12),
          legend.title = element_text(size = 15))+
    scale_fill_manual(values=c('palegreen1', 'plum1','skyblue1','dodgerblue4','black','grey40','grey80'),
                      labels = c("Shallow1","Shallow2","Deep1","Deep2"))+
    guides(fill=guide_legend(title = "Lineage", override.aes = list(shape = 21)))
legend_pca_ss_1


save(mc_bothtrans_pc12,ss_bothtrans_pc12,legend_pca_mc_1,legend_pca_mc_2,legend_pca_ss_1,file = "popgen/mc_ss_bothtrans.rdata")


#----make violin plots:------
library(tidyverse)
library(ggplot2)

spp <- 'ssid' #choose species

#import metadata file:
if(spp == 'mcav'){
  load("popgen/mcav_both_transects/mc_metadata_2ball.rdata")
}else if(spp == 'ssid'){
  load("popgen/ssid_both_transects/ss_metadata_2ball.rdata")
}

#make depth column
if(spp == 'mcav'){
  mc_meta$depth <- ifelse(mc_meta$habitat == "nearshore", 3,
                        ifelse(mc_meta$transect == "1" & mc_meta$habitat == "offshore", 7,
                               ifelse(mc_meta$transect == "2" & mc_meta$habitat == "offshore", 3,
                                      ifelse(mc_meta$transect == "1" & mc_meta$habitat == "deep", 20, 31))))
}else if(spp == 'ssid'){
  ss_meta$depth <- ifelse(ss_meta$habitat == "nearshore", 3,
                          ifelse(ss_meta$transect == "1" & ss_meta$habitat == "offshore", 7,
                                 ifelse(ss_meta$transect == "2" & ss_meta$habitat == "offshore", 3,
                                        ifelse(ss_meta$transect == "1" & ss_meta$habitat == "deep", 20, 31))))
}

#save file:
if(spp == 'mcav'){
  write.table(mc_meta, file = "mcav_files/mc_enviro.meta", col.names = T, row.names = F, quote = F)
}else if(spp == 'ssid'){
  write.table(ss_meta, file = "ssid_files/ss_enviro.meta", col.names = T, row.names = F, quote = F)
}

#rename file to env:
if(spp == 'mcav'){
  env <- mc_meta
}else if(spp == 'ssid'){
  env <- ss_meta
}

# subset to remove hybrids and juvs and put clusters in order
if(spp == 'mcav'){
  env <- subset(env, !(admix %in% c("NO1.2","OD2.3","OD2.4")))
  env2 <- env %>%
    mutate( admix=factor(admix,levels=c('N1','O2','D3','D4')) ) %>%
    subset(!age %in% c("juvenile"))
}else if(spp == 'ssid'){
  env <- subset(env, !(admix %in% c("S1.D1","S1.S2","S2.D1")))
  env2 <- env %>%
    mutate( admix=factor(admix,levels=c('S1','S2','D1','D2')) ) %>%
    subset(!age %in% c("juvenile"))
}

# Set colors for plotting
cols=c('palegreen1','plum1','skyblue','dodgerblue4')

# Violin plot (rename based on species)
if(spp=='mcav'){
  mc_violin <- ggplot(env2, aes(x = admix, y = depth, fill = admix)) +
    geom_violin(scale = "width", show.legend = F) +
    scale_fill_manual(values = cols) +
    xlab('Lineage') +
    ylab('Depth')+
    theme_classic()+
    theme(axis.title = element_text(size =16),axis.text = element_text(size = 14), axis.ticks.x = element_blank())+
    scale_x_discrete(labels = NULL)
  mc_violin
}else if(spp=='ssid'){
  ss_violin <- ggplot(env2, aes(x = admix, y = depth, fill = admix)) +
    geom_violin(scale = "width", show.legend = F) +
    scale_fill_manual(values = cols) +
    xlab('Lineage') +
    ylab('Depth')+
    theme_classic()+
    theme(axis.title = element_text(size =16),axis.text = element_text(size = 14), axis.ticks.x = element_blank())+
    scale_x_discrete(labels = NULL)
  ss_violin
}

#save ggplot data
if(spp == 'mcav'){
  save(mc_violin,file = "popgen/mcav_violinplot.rdata")
}else if(spp == 'ssid'){
  save(ss_violin, file = "popgen/ssid_violinplot.rdata")
}

# 2bRAD_set2 --------------------------------------------------------------
library(vegan)
library(ggplot2)
library(stringr)
library(dplyr)

spp <- 'ssid' #choose species
npops <- 4 #choose number of populations

dir <- paste0('popgen/',spp,'_set2/') # change this to where your scp'd files are
bams <- read.table(paste0(dir, "bams"))[,1] # list of bam files
if(spp=='ssid'){
  prefix <- 'ss'
  bams <- gsub(".nosymbio.fastq.bam","",bams)
} else if (spp=='mcav'){
  prefix <- 'mc'
  bams <- gsub("_s2.trim.bt2.bam","",bams)
}

ibsMat <- as.matrix(read.table(paste0(dir,prefix,"02b.ibsMat")))
dimnames(ibsMat) <- list(bams,bams)

hc <- hclust(as.dist(ibsMat),"ave")
plot(hc,cex=0.7) # clustering of samples by IBS (great to detect clones or closely related individuals)

if(spp=='mcav'){
  clones <- c("MJ6D", "MA6N")
}else if(spp=='ssid'){
  clones <- c("SA16D_s2", "SJ17N-r_s2", "SA11N_s2", "SA10N_s2", "SJ2Ob_s2", "SJ8Ob_s2", "SJ1N_s2", 
               "SJ18N_s2","SJ12N_s2", "SJ13N_s2","SJ19N_s2", "SJ11N_s2")
}

#remove clones
goods <- which(!(bams %in% clones))
clonn <- which((bams %in% clones))
ibsMat <- ibsMat[goods,goods]
bams <- bams[goods]

if(spp=='ssid'){
  bamlist <- data.frame(paste(bams,".nosymbio.fastq.bam",sep=""))
} else if (spp=='mcav'){
  bamlist <- data.frame(paste(bams,"_s2.trim.bt2.bam",sep=""))
}
write.table(bamlist,file=paste0(dir,"bams_noclones"),row.names=F,col.names=F,quote=F) # use this bam list to rerun ANGSD for ADMIXTURE

#---------------
#after rerunning ANGSD for ADMIXTURE
#clear environment before redoing with new IBS matrix without clones
spp <- 'mcav' #choose species
npops <- 4 #choose number of populations

# set 2 2brad without clones ----------------------------------------------

dir <- paste0('popgen/',spp,'_set2/') # change this to where your scp'd files are
bams <- read.table(paste0(dir, "bams_noclones"))[,1] # list of bam files
if(spp=='ssid'){
  prefix <- 'ss'
  bams <- gsub("_s2.nosymbio.fastq.bam","",bams)
} else if (spp=='mcav'){
  prefix <- 'mc'
  bams <- gsub("_s2.trim.bt2.bam","",bams)
}

ibsMat <- as.matrix(read.table(paste0(dir,spp,"2b.ibsMat")))
dimnames(ibsMat) <- list(bams,bams)

hc <- hclust(as.dist(ibsMat),"ave")
plot(hc,cex=0.7) # clustering of samples by IBS (great to detect clones or closely related individuals)

save(ibsMat,bams, file = paste0(dir,spp,"_2b_set2_ibs_bams.rdata"))

#==================================
# Make metadata files for set 2 ----------------------------------------------------------
spp <- 'ssid' #choose species
npops <- 4 #choose number of populations

if(spp=='ssid'){
  prefix <- 'ss'
} else if (spp=='mcav'){
  prefix <- 'mc'
}

#load ibs and bams:
load(paste0("popgen/",spp,"_set2/",spp,"_2b_set2_ibs_bams.rdata"))

adults=grep("A",bams)
juveniles=grep("J",bams)
near=grep("N",bams)
off=grep("O",bams)
deep=grep("D",bams)

# factors:
site=rep("deep",length(bams))
site[off]="offshore"
site[near]="nearshore"
aj=rep("adult",length(bams))
aj[juveniles]="juvenile"

load(paste0('popgen/',spp,"_set2/",prefix,"_k4.qopt_clusters.RData"))

if(spp=='mcav'){
  mc_meta <- data.frame(id = bams, age = aj, habitat = site, admix = cluster.admix)
  mc_meta <- mc_meta %>%
    mutate(admix = recode(admix, "1"="N1", "2"="D3", "3"="O2", "4"="OD2.3", "1.3"="NO1.2", 
                          "1.4"="D4","2.3"="OD2.4", "3.4"="OD2.4"))
  mc_meta["MA12O","admix"] <- "OD2.4"
  mc_meta["MJ7D","admix"] <- "D4"
}else if(spp=='ssid'){
  ss_meta <- data.frame(id = bams, age = aj, habitat = site, admix = cluster.admix)
  ss_meta <- ss_meta %>%
    mutate(admix = recode(admix, "1"="D1","1.2"="S2","2"="S1","2.3"="S1","3"="S1"))
  ss_meta["SA3D","admix"] <- "D2"
  ss_meta["SA7D","admix"] <- "D2"
  ss_meta["SJ5O","admix"] <- "S2.D1"
}

#save file
if(spp=='mcav'){
  write.table(mc_meta,file="mcav_files/mcav.metadata.2b_set2",row.names=F,col.names=F,quote=F)
  save(mc_meta,file = "popgen/mcav_set2/mc_metadata_2b_set2.rdata")
}else if(spp=='ssid'){
  write.table(ss_meta,file="ssid_files/ssid.metadata.2b_set2",row.names=F,col.names=F,quote=F)
  save(ss_meta,file = "popgen/ssid_set2/ss_metadata_2b_set2.rdata")
}
