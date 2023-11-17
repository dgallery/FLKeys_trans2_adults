### GO-MWU on ssid adults WGCNA
setwd("ssid_go/")
#################
goAnnotations="ssid_go.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
source("gomwu.functions.R")
MFgoDivision="MF" 
BPgoDivision="BP"
CCgoDivision="CC"

#turquoise module
turquoiseinput="turquoise.csv" 

#turquoise MF stats/results:
gomwuStats(turquoiseinput, goDatabase, goAnnotations, MFgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
turquoiseMFresults=gomwuPlot(turquoiseinput,goAnnotations,MFgoDivision,
                         absValue = 0.001,
                         level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                         level2=0.05, # FDR cutoff to print in regular (not italic) font.
                         level3=0.01, # FDR cutoff to print in large bold font.
                         txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                         treeHeight=0.5, # height of the hierarchical clustering tree
)
turquoiseMFresults

#turquoise BP stats/results:
gomwuStats(turquoiseinput, goDatabase, goAnnotations, BPgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered (150, 100 and 50 for BP)
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
turquoiseBPresults=gomwuPlot(turquoiseinput,goAnnotations,BPgoDivision,
                         absValue = 0.001,
                         level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                         level2=0.05, # FDR cutoff to print in regular (not italic) font.
                         level3=0.01, # FDR cutoff to print in large bold font.
                         txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                         treeHeight=0.5, # height of the hierarchical clustering tree
)
turquoiseBPresults  

#turquoise CC stats/results:
gomwuStats(turquoiseinput, goDatabase, goAnnotations, CCgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
turquoiseCCresults=gomwuPlot(turquoiseinput,goAnnotations,CCgoDivision,
                         absValue = 0.001,
                         level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                         level2=0.05, # FDR cutoff to print in regular (not italic) font.
                         level3=0.01, # FDR cutoff to print in large bold font.
                         txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                         treeHeight=0.5, # height of the hierarchical clustering tree
)
turquoiseCCresults  

#green module
greeninput="green.csv" 

#green MF stats/results:
gomwuStats(greeninput, goDatabase, goAnnotations, MFgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
greenMFresults=gomwuPlot(greeninput,goAnnotations,MFgoDivision,
                        absValue = 0.001,
                        level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
greenMFresults

#green BP stats/results:
gomwuStats(greeninput, goDatabase, goAnnotations, BPgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered (150, 100 and 50 for BP)
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
greenBPresults=gomwuPlot(greeninput,goAnnotations,BPgoDivision,
                        absValue = 0.001,
                        level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
greenBPresults  

#green CC stats/results:
gomwuStats(greeninput, goDatabase, goAnnotations, CCgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
greenCCresults=gomwuPlot(greeninput,goAnnotations,CCgoDivision,
                        absValue = 0.001,
                        level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
greenCCresults  

#lightyellow module
lightyellowinput="lightyellow.csv" 

#lightyellow MF stats/results:
gomwuStats(lightyellowinput, goDatabase, goAnnotations, MFgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
lightyellowMFresults=gomwuPlot(lightyellowinput,goAnnotations,MFgoDivision,
                             absValue = 0.001,
                             level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                             level2=0.05, # FDR cutoff to print in regular (not italic) font.
                             level3=0.01, # FDR cutoff to print in large bold font.
                             txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                             treeHeight=0.5, # height of the hierarchical clustering tree
)
lightyellowMFresults

#lightyellow BP stats/results:
gomwuStats(lightyellowinput, goDatabase, goAnnotations, BPgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered (150, 100 and 50 for BP)
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
lightyellowBPresults=gomwuPlot(lightyellowinput,goAnnotations,BPgoDivision,
                             absValue = 0.001,
                             level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                             level2=0.05, # FDR cutoff to print in regular (not italic) font.
                             level3=0.01, # FDR cutoff to print in large bold font.
                             txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                             treeHeight=0.5, # height of the hierarchical clustering tree
)
lightyellowBPresults  

#lightyellow CC stats/results:
gomwuStats(lightyellowinput, goDatabase, goAnnotations, CCgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
lightyellowCCresults=gomwuPlot(lightyellowinput,goAnnotations,CCgoDivision,
                             absValue = 0.001,
                             level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                             level2=0.05, # FDR cutoff to print in regular (not italic) font.
                             level3=0.01, # FDR cutoff to print in large bold font.
                             txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                             treeHeight=0.5, # height of the hierarchical clustering tree
)
lightyellowCCresults

#royalblue module
royalblueinput="royalblue.csv" 

#royalblue MF stats/results:
gomwuStats(royalblueinput, goDatabase, goAnnotations, MFgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
royalblueMFresults=gomwuPlot(royalblueinput,goAnnotations,MFgoDivision,
                             absValue = 0.001,
                             level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                             level2=0.05, # FDR cutoff to print in regular (not italic) font.
                             level3=0.01, # FDR cutoff to print in large bold font.
                             txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                             treeHeight=0.5, # height of the hierarchical clustering tree
)
royalblueMFresults

#royalblue BP stats/results:
gomwuStats(royalblueinput, goDatabase, goAnnotations, BPgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered (150, 100 and 50 for BP)
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
royalblueBPresults=gomwuPlot(royalblueinput,goAnnotations,BPgoDivision,
                             absValue = 0.001,
                             level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                             level2=0.05, # FDR cutoff to print in regular (not italic) font.
                             level3=0.01, # FDR cutoff to print in large bold font.
                             txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                             treeHeight=0.5, # height of the hierarchical clustering tree
)
royalblueBPresults  

#royalblue CC stats/results:
gomwuStats(royalblueinput, goDatabase, goAnnotations, CCgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
royalblueCCresults=gomwuPlot(royalblueinput,goAnnotations,CCgoDivision,
                             absValue = 0.001,
                             level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                             level2=0.05, # FDR cutoff to print in regular (not italic) font.
                             level3=0.01, # FDR cutoff to print in large bold font.
                             txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                             treeHeight=0.5, # height of the hierarchical clustering tree
)
royalblueCCresults

#greenyellow module
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
           smallest=5,   # a GO category should contain at least this many genes to be considered (150, 100 and 50 for BP)
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

#grey60 module
grey60input="grey60.csv" 

#grey60 MF stats/results:
gomwuStats(grey60input, goDatabase, goAnnotations, MFgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
grey60MFresults=gomwuPlot(grey60input,goAnnotations,MFgoDivision,
                        absValue = 0.001,
                        level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
grey60MFresults

#grey60 BP stats/results:
gomwuStats(grey60input, goDatabase, goAnnotations, BPgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered (150, 100 and 50 for BP)
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
grey60BPresults=gomwuPlot(grey60input,goAnnotations,BPgoDivision,
                        absValue = 0.001,
                        level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
grey60BPresults  

#grey60 CC stats/results:
gomwuStats(grey60input, goDatabase, goAnnotations, CCgoDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
)
grey60CCresults=gomwuPlot(grey60input,goAnnotations,CCgoDivision,
                        absValue = 0.001,
                        level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                        level2=0.05, # FDR cutoff to print in regular (not italic) font.
                        level3=0.01, # FDR cutoff to print in large bold font.
                        txtsize=1.5,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                        treeHeight=0.5, # height of the hierarchical clustering tree
)
grey60CCresults

save(turquoiseMFresults, turquoiseBPresults, turquoiseCCresults, greenMFresults, greenBPresults, greenCCresults,
     royalblueMFresults, royalblueBPresults, royalblueCCresults, greenyellowMFresults, greenyellowBPresults, 
     grey60MFresults, grey60BPresults, grey60CCresults, file = "ssid_go_results.rdata")

##### PLOTS #######
setwd("ssid_go/")
load("ssid_go_results.rdata")

GO_terms_list <- list(turquoiseMFresults, turquoiseBPresults, turquoiseCCresults, greenMFresults, greenBPresults, greenCCresults,
                      royalblueMFresults, royalblueBPresults, royalblueCCresults, greenyellowMFresults, greenyellowBPresults, 
                      grey60MFresults, grey60BPresults, grey60CCresults)

names(GO_terms_list) <- c("turquoiseMFresults", "turquoiseBPresults", "turquoiseCCresults", "greenMFresults", "greenBPresults", "greenCCresults",
                          "royalblueMFresults", "royalblueBPresults", "royalblueCCresults", "greenyellowMFresults", "greenyellowBPresults", 
                          "grey60MFresults", "grey60BPresults", "grey60CCresults")

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
       tip.color = GO_color_list[[i]][[2]]$newcolor,
       main = names(GO_color_list[i]),
       cex = GO_color_list[[i]][[1]]$text,
       font = GO_color_list[[i]][[1]]$face,
       label.offset = 0.2)
}

#plot individual results one at a time
par(mar=c(1,0.5,3,0.5))
module <- "turquoise"
division <- "BP"
plot(as.phylo(GO_color_list[[paste0(module,division,"results")]][[2]]),
     tip.color = GO_color_list[[paste0(module,division,"results")]][[2]]$color,
     main = names(GO_color_list[paste0(module,division,"results")]),
     cex = GO_color_list[[paste0(module,division,"results")]][[1]]$text,
     font = GO_color_list[[paste0(module,division,"results")]][[1]]$face,
     use.edge.length = F)
legend(x = "topright", text.font = c(2,3,1), text.col = "black", legend=c(expression(bold("p"<="0.01")),expression(italic("p"<="0.05")),expression("p"<="0.1")))

