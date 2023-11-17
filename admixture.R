#~~~~~~
# Code for assigning admixture (modified from JP's plot_admixture.R file)
#~~~~~~

#---------------
spp <- 'mcav' #choose species
npops <- 4 #choose number of populations
#---------------

# assembling the input table

#dir=paste0('popgen/',spp,"_both_transects/") # path to input files - both sets
dir=paste0('popgen/',spp,"_set2/") # path to input files - set 2

if(spp=='ssid'){
  prefix <- 'ss'
  inName <- paste0('ss_k', npops, '.qopt')
} else if (spp=='mcav'){
  prefix <- 'mc'
  inName <- paste0('mc_k', npops, '.qopt')
}

#pops="inds2pops_all" # both sets - 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
pops = "inds2pops" #set 2 only

#------------

npops=as.numeric(sub("\\D+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(dir,inName,sep=""),header=F)
i2p=read.table(paste(dir,pops,sep=""),header=F)
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
rownames(tbl) <- tbl$ind <- sub("(.*?)\\..*$", "\\1", tbl$ind)
head(tbl,20) # this is how the resulting dataset must look

# recording cluster affiliations
cluster.admix=apply(tbl[,1:npops],1,function(x) {return(paste(which(x>0.25),collapse=".")) })
save(cluster.admix,file=paste0(dir,inName,"_clusters.RData"))
