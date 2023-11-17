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
axes2plot=c(1,3)
axes=1.3

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
} else if (spp == 'ssid'){
  admix.spider$admix <- factor(admix.spider$admix, levels = c("S1","S2","D1","D2")) #no hybrids in adult ssids
}

#ggplot
ggplot(data=admix.spider)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    labs(x = paste0('MDS', axes2plot[1], ' (', round(eigvals[2,axes2plot[1]]*100,1), '%)'), y = paste0('MDS', axes2plot[2], ' (', round(eigvals[2,axes2plot[2]]*100,1), '%)'))+
    geom_segment(aes_string(x='mean.x', y='mean.y', xend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], yend=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], alpha='transect'), col='grey70')+
    geom_point(aes_string(x=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[1]], y=colnames(scores(pp0_ad, display = 'sites', 1:4))[axes2plot[2]], shape='transect', fill='admix', col='transect', alpha='transect'), size = 5)+
    scale_color_manual(values = c('black','black'))+
    scale_alpha_manual(values = c(1,0.8)) +
    scale_fill_manual(values = c('palegreen1', 'plum1','skyblue1','dodgerblue4','black','grey40','grey80'))+
    scale_shape_manual(aes(transect),values = c(21,24))+
    #coord_fixed(xlim = c(-0.15,0.3))+ #ssid
    #coord_fixed(xlim = c(-0.22,0.1))+  #mcav pc1,2
    theme(legend.position = "none", axis.title = element_text(size =16),axis.text = element_text(size = 14))
  #guides(fill=guide_legend(override.aes = list(colour=c('palegreen1', 'plum1','skyblue1','dodgerblue4'))), shape = "none", alpha = "none", col = "none")

