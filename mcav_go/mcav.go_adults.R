#Go_MWU for mcav WGCNA modules 0.5 eigengene merge from 9-8-2021
setwd("mcav_go/")

##----------Normal GO-MWU for within species comparisons ####
goAnnotations="mcav_go.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
source("gomwu.functions.R")
MFgoDivision="MF" 
BPgoDivision="BP"
CCgoDivision="CC"

#lightcyan module
lightcyaninput="lightcyan.csv" 

#lightcyan MF stats/results:
gomwuStats(lightcyaninput, goDatabase, goAnnotations, MFgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
lightcyanMFresults=gomwuPlot(lightcyaninput,goAnnotations,MFgoDivision,
                         absValue = 0.001,
                         level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                         level2=0.05, # FDR cutoff to print in regular (not italic) font.
                         level3=0.01, # FDR cutoff to print in large bold font.
                         txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                         treeHeight=0.5, # height of the hierarchical clustering tree
)
lightcyanMFresults

#lightcyan BP stats/results:
gomwuStats(lightcyaninput, goDatabase, goAnnotations, BPgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered (150, 100 and 50 for BP)
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
lightcyanBPresults=gomwuPlot(lightcyaninput,goAnnotations,BPgoDivision,
                         absValue = 0.001,
                         level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                         level2=0.05, # FDR cutoff to print in regular (not italic) font.
                         level3=0.01, # FDR cutoff to print in large bold font.
                         txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                         treeHeight=0.5, # height of the hierarchical clustering tree
)
lightcyanBPresults  

#lightcyan CC stats/results:
gomwuStats(lightcyaninput, goDatabase, goAnnotations, CCgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
lightcyanCCresults=gomwuPlot(lightcyaninput,goAnnotations,CCgoDivision,
                         absValue = 0.001,
                         level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                         level2=0.05, # FDR cutoff to print in regular (not italic) font.
                         level3=0.01, # FDR cutoff to print in large bold font.
                         txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                         treeHeight=0.5, # height of the hierarchical clustering tree
)
lightcyanCCresults  

#lightgreen module
lightgreeninput="lightgreen.csv" 

#lightgreen MF stats/results:
gomwuStats(lightgreeninput, goDatabase, goAnnotations, MFgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
lightgreenMFresults=gomwuPlot(lightgreeninput,goAnnotations,MFgoDivision,
                             absValue = 0.001,
                             level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                             level2=0.05, # FDR cutoff to print in regular (not italic) font.
                             level3=0.01, # FDR cutoff to print in large bold font.
                             txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                             treeHeight=0.5, # height of the hierarchical clustering tree
)
lightgreenMFresults  

#lightgreen BP stats/results:
gomwuStats(lightgreeninput, goDatabase, goAnnotations, BPgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
lightgreenBPresults=gomwuPlot(lightgreeninput,goAnnotations,BPgoDivision,
                             absValue = 0.001,
                             level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                             level2=0.05, # FDR cutoff to print in regular (not italic) font.
                             level3=0.01, # FDR cutoff to print in large bold font.
                             txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                             treeHeight=0.5, # height of the hierarchical clustering tree
)
lightgreenBPresults 

#lightgreen CC stats/results:
gomwuStats(lightgreeninput, goDatabase, goAnnotations, CCgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
lightgreenCCresults=gomwuPlot(lightgreeninput,goAnnotations,CCgoDivision,
                             absValue = 0.001,
                             level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                             level2=0.05, # FDR cutoff to print in regular (not italic) font.
                             level3=0.01, # FDR cutoff to print in large bold font.
                             txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                             treeHeight=0.5, # height of the hierarchical clustering tree
)
lightgreenCCresults  

#darkgrey module
darkgreyinput="darkgrey.csv" 

#darkgrey MF stats/results:
gomwuStats(darkgreyinput, goDatabase, goAnnotations, MFgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be consideturquoise if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be consideturquoise
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
darkgreyMFresults=gomwuPlot(darkgreyinput,goAnnotations,MFgoDivision,
                           absValue = 0.001,
                           level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                           level2=0.05, # FDR cutoff to print in regular (not italic) font.
                           level3=0.01, # FDR cutoff to print in large bold font.
                           txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                           treeHeight=0.5, # height of the hierarchical clustering tree
)
darkgreyMFresults  

#darkgrey BP stats/results:
gomwuStats(darkgreyinput, goDatabase, goAnnotations, BPgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be consideturquoise if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be consideturquoise
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
darkgreyBPresults=gomwuPlot(darkgreyinput,goAnnotations,BPgoDivision,
                           absValue = 0.001,
                           level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                           level2=0.05, # FDR cutoff to print in regular (not italic) font.
                           level3=0.01, # FDR cutoff to print in large bold font.
                           txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                           treeHeight=0.5, # height of the hierarchical clustering tree
)
darkgreyBPresults  

#darkgrey CC stats/results:
gomwuStats(darkgreyinput, goDatabase, goAnnotations, CCgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be consideturquoise if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be consideturquoise
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
darkgreyCCresults=gomwuPlot(darkgreyinput,goAnnotations,CCgoDivision,
                           absValue = 0.001,
                           level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                           level2=0.05, # FDR cutoff to print in regular (not italic) font.
                           level3=0.01, # FDR cutoff to print in large bold font.
                           txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                           treeHeight=0.5, # height of the hierarchical clustering tree
)
darkgreyCCresults  

#darkturquoise module
darkturquoiseinput="darkturquoise.csv" 

#darkturquoise MF stats/results:
gomwuStats(darkturquoiseinput, goDatabase, goAnnotations, MFgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
darkturquoiseMFresults=gomwuPlot(darkturquoiseinput,goAnnotations,MFgoDivision,
                        absValue = 0.001,
                        level1=0.5, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
darkturquoiseMFresults  

#darkturquoise BP stats/results:
gomwuStats(darkturquoiseinput, goDatabase, goAnnotations, BPgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered (still 37 after bumping up to 0.05 - bump this to 100)
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
darkturquoiseBPresults=gomwuPlot(darkturquoiseinput,goAnnotations,BPgoDivision,
                        absValue = 0.001,
                        level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
darkturquoiseBPresults  

#darkturquoise CC stats/results:
gomwuStats(darkturquoiseinput, goDatabase, goAnnotations, CCgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
darkturquoiseCCresults=gomwuPlot(darkturquoiseinput,goAnnotations,CCgoDivision,
                        absValue = 0.001,
                        level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
darkturquoiseCCresults  

#midnightblue module
greenyellowinput="greenyellow.csv" 

#greenyellow MF stats/results:
gomwuStats(greenyellowinput, goDatabase, goAnnotations, MFgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
greenyellowMFresults=gomwuPlot(greenyellowinput,goAnnotations,MFgoDivision,
                        absValue = 0.001,
                        level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
greenyellowMFresults  

#greenyellow BP stats/results:
gomwuStats(greenyellowinput, goDatabase, goAnnotations, BPgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
greenyellowBPresults=gomwuPlot(greenyellowinput,goAnnotations,BPgoDivision,
                        absValue = 0.001,
                        level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
greenyellowBPresults  

#greenyellow CC stats/results:
gomwuStats(greenyellowinput, goDatabase, goAnnotations, CCgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
greenyellowCCresults=gomwuPlot(greenyellowinput,goAnnotations,CCgoDivision,
                        absValue = 0.001,
                        level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
greenyellowCCresults  

#midnightblue module
midnightblueinput="midnightblue.csv" 

#midnightblue MF stats/results:
gomwuStats(midnightblueinput, goDatabase, goAnnotations, MFgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
midnightblueMFresults=gomwuPlot(midnightblueinput,goAnnotations,MFgoDivision,
                               absValue = 0.001,
                               level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                               level2=0.05, # FDR cutoff to print in regular (not italic) font.
                               level3=0.01, # FDR cutoff to print in large bold font.
                               txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                               treeHeight=0.5, # height of the hierarchical clustering tree
)
midnightblueMFresults  

#midnightblue BP stats/results:
gomwuStats(midnightblueinput, goDatabase, goAnnotations, BPgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
midnightblueBPresults=gomwuPlot(midnightblueinput,goAnnotations,BPgoDivision,
                               absValue = 0.001,
                               level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                               level2=0.05, # FDR cutoff to print in regular (not italic) font.
                               level3=0.01, # FDR cutoff to print in large bold font.
                               txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                               treeHeight=0.5, # height of the hierarchical clustering tree
)
midnightblueBPresults  

#midnightblue CC stats/results:
gomwuStats(midnightblueinput, goDatabase, goAnnotations, CCgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
midnightblueCCresults=gomwuPlot(midnightblueinput,goAnnotations,CCgoDivision,
                               absValue = 0.001,
                               level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                               level2=0.05, # FDR cutoff to print in regular (not italic) font.
                               level3=0.01, # FDR cutoff to print in large bold font.
                               txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                               treeHeight=0.5, # height of the hierarchical clustering tree
)
midnightblueCCresults  

save(lightcyanMFresults,lightcyanBPresults,lightcyanCCresults, lightgreenMFresults,lightgreenBPresults,lightgreenCCresults, darkgreyBPresults, 
     darkgreyCCresults, greenyellowMFresults, greenyellowBPresults, greenyellowCCresults, file = "mcav_go_results.rdata")

##### PLOTS #######
setwd("mcav_go/")
load("mcav_go_results.rdata")
library(ape)

GO_terms_list <- list(lightcyanMFresults,lightcyanBPresults,lightcyanCCresults, lightgreenMFresults,lightgreenBPresults,lightgreenCCresults, darkgreyBPresults, 
                      darkgreyCCresults, greenyellowMFresults, greenyellowBPresults, greenyellowCCresults)

names(GO_terms_list) <- c("lightcyanMFresults","lightcyanBPresults","lightcyanCCresults", "lightgreenMFresults","lightgreenBPresults","lightgreenCCresults", "darkgreyBPresults", 
                          "darkgreyCCresults", "greenyellowMFresults", "greenyellowBPresults", "greenyellowCCresults")

GO_color_list <- GO_terms_list

pval <- matrix(nrow=0,ncol=2)

for (i in seq_along(GO_color_list)) {
  GO_color_list[[i]][[2]]$newcolor <- 
    ifelse(grepl("translat|ribosom|rRNA|ribonu",GO_color_list[[i]][[2]]$labels,ignore.case = T),"brown",
           ifelse(grepl("actin|bind|chromo|condens|cytoplas",GO_color_list[[i]][[2]]$labels,ignore.case = T),"blue",
                  ifelse(grepl("metab|proteo|peptidase|secretory|catab",GO_color_list[[i]][[2]]$labels,ignore.case = T),"green4","black")))
  pval <- cbind(row.names(GO_color_list[[i]][[1]]),GO_color_list[[i]][[1]]$pval)
  reorder_idx <- match(GO_color_list[[i]][[2]]$labels,pval)
  pval2 <- pval[,2]
  pval2 <- pval2[reorder_idx]
  GO_color_list[[i]][[1]] <- GO_color_list[[i]][[1]][match(pval2,GO_color_list[[i]][[1]]$pval),]
  #GO_color_list[[i]][[2]]$pval <- pval2
  GO_color_list[[i]][[1]]$text <-
    ifelse(GO_color_list[[i]][[1]]$pval > 0.05, 0.8,
           ifelse(GO_color_list[[i]][[1]]$pval <=0.01, 1.2, 1))
  GO_color_list[[i]][[1]]$face <-
    ifelse(GO_color_list[[i]][[1]]$pval > 0.05, 1,
           ifelse(GO_color_list[[i]][[1]]$pval <=0.01, 2, 3))
}

for (i in seq_along(GO_color_list)) {
  plot(as.phylo(GO_color_list[[i]][[2]]),
       tip.color = GO_color_list[[i]][[2]]$color,
       main = names(GO_color_list[i]),
       cex = GO_color_list[[i]][[1]]$text,
       font = GO_color_list[[i]][[1]]$face,
       label.offset = 0.2)
}

#plot individual results one at a time
par(mar=c(1,0.5,3,0.5))
module <- "darkgrey"
division <- "BP"
plot(as.phylo(GO_color_list[[paste0(module,division,"results")]][[2]]),
     tip.color = GO_color_list[[paste0(module,division,"results")]][[2]]$color,
     main = names(GO_color_list[paste0(module,division,"results")]),
     cex = GO_color_list[[paste0(module,division,"results")]][[1]]$text,
     font = GO_color_list[[paste0(module,division,"results")]][[1]]$face,
     use.edge.length = F)

# text(40,9,substitute(bold("p"<= "0.01")),font=GO_color_list[[paste0(module,division,"results")]][[1]]$face,
#      cex=GO_color_list[[paste0(module,division,"results")]][[1]]$text)
# text(40,8,substitute(bold("p"<= "0.05")),font=GO_color_list[[paste0(module,division,"results")]][[1]]$face,
#      cex=GO_color_list[[paste0(module,division,"results")]][[1]]$text)
# text(40,7,substitute(bold("p"<= "0.05")),font=GO_color_list[[paste0(module,division,"results")]][[1]]$face,
#      cex=GO_color_list[[paste0(module,division,"results")]][[1]]$text)

legend(x = "topright", text.font = c(2,3,1), text.col = "black", legend=c(expression(bold("p"<="0.01")),expression(italic("p"<="0.05")),expression("p"<="0.1")))
