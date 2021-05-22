# Run edgeR on PhenoSOM results
# Author:  Erin Simonds
# Version: 1.210522
#
#
# Overall workflow of this script:   
#   Read in the results from a PhenoSOM run
#	Read in a setup CSV file defining conditions A and  B
#   Compare cluster abundance in condition B relative to A
#   Note: The baseline condition (Condition A) will be whichever condition comes first in the alphabet. It can be helpful to prepend the control condition with a period (e.g. ".Vehicle") to make sure it is always used as the baseline.
#
# Input:
#   1. A subfolder containing the results of a PhenoSOM run.
#   2. A CSV file containing a key of the conditions. The filename must end with "setup.csv"
#
# Output:
#   A subfolder containing edgeR results as heatmaps, boxplots, and a PCA plot of metaclusters


rm(list = ls())

library(edgeR)
library(scran)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(tidyr)

################################################# User-defined parameters ##################################################

# Define where to find the PhenoSOM clustering output
clusteringdir <- "PhenoSOM_Step2_output"

# Define a general output folder for all edgeR analyses
outdir <- "PhenoSOM_Step3_output"

# Define a name for this edgeR analysis (this will be appended to all output filenames)
analysisname <- "What conditions are you comparing?"

############################################################################################################################

# Set working directory to the location of the script, which has been Sourced from within Rstudio
setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path)) # Set working directory to location of currently open Rstudio script
startdir <- getwd()
set.seed(42)

# Create output folders
outdir2 <- file.path(outdir, analysisname)
dir.create(outdir, showWarnings = F)
dir.create(outdir2, showWarnings = F)

# Check if edgeR setup CSV file exists
if(length(dir(path=".", pattern=".*setup\\.csv")) == 0){
  stop("[Error] Can't find CSV file ")
}

# Set the conditions we are comparing (i.e. A vs. B)
setup_table <- read.csv(dir(path=".", pattern=".*setup\\.csv")[1], row.names=NULL)

# Make sure setup_table contains the required columns
if(all(c("FCSGUID", "Condition") %in% colnames(setup_table))){
  print("[OK] edgeR setup CSV file loaded")
} else {
  stop("[Error] edgeR setup CSV file must contain the columns 'FCSGUID' and 'Condition'")
}

# Check if clustering folder exists
if(length(dir(clusteringdir)) == 0){
  stop("[Error] Can't find clustering output folder")
}

# Load old clustering data
if(length(dir(clusteringdir, pattern=".*\\.Rdata")) != 1){
  stop("[Error] Can't find Rdata file from clustering output")
} else {
  oldprojectname <- gsub(pattern="\\ results\\.Rdata", x=dir(clusteringdir, pattern=".*\\.Rdata"), replacement = "")  
}
key.channels <- read.csv(file=file.path(clusteringdir, paste0(oldprojectname, " channel and marker names used.csv")))
plottingChannels <- as.character(key.channels$MarkerChannel[as.logical(key.channels$IsPlottingChannel)])
clusteringChannels <- as.character(key.channels$MarkerChannel[as.logical(key.channels$IsClusteringChannel)])
key.fcsids <- read.csv(file=file.path(clusteringdir, paste0(oldprojectname, " sample IDs and sizes in input dataset.csv", row.names=NULL)))
clustercentroids <- read.csv(file=file.path(clusteringdir, "Maps and plots of PhenoGraph clusters", paste0(oldprojectname, " Rphenograph metacluster medians - all markers.csv", row.names=NULL)))

## Sanity checks
if(all(length(key.fcsids$FCSGUID) == length(setup_table$FCSGUID))){
  print("[OK] Number of FCS files in conditions setup table matches key")
} else {
  stop("[Error] Number of FCS files in conditions setup table does not match key")
}

if(all(sort(unique(key.fcsids$FCSGUID)) == sort(unique(setup_table$FCSGUID)))){
  print("[OK] FCS filenames in conditions setup table match key")
} else {
  stop("[Error] FCS filenames in conditions setup table do not match key")
}

# extract the conditions we are comparing, ordering the vector by the FCSGUID in key.fcsid
conditions <- setup_table$Condition[match(key.fcsids$FCSGUID, setup_table$FCSGUID)]
has_condition_mask <- (!is.na(conditions))
print(table(conditions))

colnames(clustercentroids)[-1] <- as.character(key.channels$MarkerChannel)
rownames(clustercentroids) <- sprintf("%02d", clustercentroids$PhenoGraphMetacluster)
clusterpercentages.orig <- t(read.csv(file=file.path(clusteringdir, "Maps and plots of PhenoGraph clusters", paste0(oldprojectname, " Rphenograph metacluster sizes (percentage) per file.csv")), row.names = NULL)[,-1])
colnames(clusterpercentages.orig) <- as.character(key.fcsids$FCSprettyname)
rownames(clusterpercentages.orig) <- sprintf("%02d", clustercentroids$PhenoGraphMetacluster)
clustercounts.orig <- t(read.csv(file=file.path(clusteringdir, "Maps and plots of PhenoGraph clusters", paste0(oldprojectname, " Rphenograph metacluster sizes (cell number) per file.csv")), row.names = NULL)[,-1])
colnames(clustercounts.orig) <- as.character(key.fcsids$FCSprettyname)
rownames(clustercounts.orig) <- sprintf("%02d", clustercentroids$PhenoGraphMetacluster)

# Remove files that are not being compared, then calculate cell counts and percentages in each cluster
clustercounts <- clustercounts.orig[,has_condition_mask]
clusterpercentages <- clusterpercentages.orig[,has_condition_mask]
sampleTotals <- colSums(clustercounts)
clusterAvgPercentages <- rowMeans(clusterpercentages)

# We can use a number of methods to test the count data for differential abundance. Here, we will use the quasi-likelihood (QL) method from the edgeR package. This allows us to model discrete count data with overdispersion due to biological variability.
y.orig <- DGEList(as.matrix(clustercounts), lib.size=sampleTotals)

# First, we do some filtering to remove irrelevant FCS files (condition = NA), and to remove low-abundance metaclusters with average counts below 5. These are mostly uninteresting as they do not provide enough evidence to reject the null hypothesis. Removing them also reduces computational work and the severity of the multiple testing correction. Lower values can also be used, but we do not recommend going below 1.
keep <- aveLogCPM(y.orig) >= aveLogCPM(5, mean(sampleTotals))
cd <- clustercounts[keep,]
y <- y.orig[keep,]

# We then apply the QL framework to estimate the dispersions, fit a generalized linear model and test for significant differences between conditions. We refer interested readers to the edgeR user's guide for more details.
design <- model.matrix(~factor(conditions))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit, coef=2)
results <- res$table
is.sig <- results$PValue <= 0.05
pcamap <- prcomp(clustercentroids[keep,clusteringChannels])
clusterlabels <- gsub(rownames(results), pattern=".*\\_", replacement="")

# Pretty plots of P-values vs fold-change (volcano)
plotdata <- results
xlims <- c(-max(abs(min(plotdata$logFC)), abs(max(plotdata$logFC))), max(abs(min(plotdata$logFC)), abs(max(plotdata$logFC))))
direction <- as.character(results$logFC > 0)
direction[plotdata$logFC > 0] <- "Increased"
direction[plotdata$logFC <= 0] <- "Decreased"
pdf(file = file.path(outdir2, paste0(analysisname," Rphenograph metaclusters edgeR volcano plot.pdf")), width=8, height=5)
p <- ggplot(plotdata, aes(logFC, -log10(PValue), label = clusterlabels))
print(p + geom_point(aes(colour = direction, size = clusterAvgPercentages[keep], alpha = is.sig)) + xlim(xlims) +
        scale_colour_manual("Direction", values = c("Increased" = "red4", "Decreased" = "blue4")) + xlab(paste0("log2(fold-change(", levels(factor(conditions))[2], " / ", levels(factor(conditions))[1], ")")) + ylab("-log10(P-value)") +
        scale_size_area("Avg. % of sample", max_size = 8) + geom_text(vjust = 0, nudge_y = max(-log10(plotdata$PValue))/20) +
        scale_alpha_discrete("P < 0.05", range = c(0.2, 1.0)) + 
        ggtitle(analysisname, subtitle = paste("Rphenograph metaclusters edgeR volcano plot")))
dev.off()

# Pretty plots of P-values in PCA space
plotdata <- cbind(results, pcamap$x[,1:2])
direction <- as.character(results$logFC > 0)
direction[plotdata$logFC > 0] <- "Increased"
direction[plotdata$logFC <= 0] <- "Decreased"
pdf(file = file.path(outdir2, paste0(analysisname," Rphenograph metaclusters edgeR PCA plot.pdf")), width=8, height=5)
p <- ggplot(plotdata, aes(PC1, PC2, label = clusterlabels))
print(p + geom_point(aes(colour = direction, size = clusterAvgPercentages[keep], alpha = is.sig)) +
        scale_colour_manual("Direction", values = c("Increased" = "red4", "Decreased" = "blue4")) +
        scale_size_area("Avg. % of sample", max_size = 8) + geom_text(vjust = 0, nudge_y = 0.1) +
        scale_alpha_discrete("P < 0.05", range = c(0.2, 1.0)) + 
        ggtitle(analysisname, subtitle = paste("Rphenograph metaclusters in PCA space"))
)
dev.off()

# Heatmap of metacluster medians, with signifant metaclusters noted
plotdata <- results
mycolorfun <- colorRampPalette(rev(brewer.pal(11,"RdBu")))
negFCsigextreme <- min(-1e7, plotdata$logFC[is.sig][plotdata$logFC[is.sig] < 0])
negFCsigthreshold <- max(-1e6, plotdata$logFC[is.sig][plotdata$logFC[is.sig] < 0])
posFCsigthreshold <- min(1e6, plotdata$logFC[is.sig][plotdata$logFC[is.sig] > 0])
posFCsigextreme <- max(1e7, plotdata$logFC[is.sig][plotdata$logFC[is.sig] > 0])
FCbreaks <- c(negFCsigextreme-1, negFCsigthreshold+.0001, posFCsigthreshold-.0001, posFCsigextreme+1)
clustercolors <- c("blue4", "gray80", "red4")[cut(plotdata$logFC,breaks=FCbreaks)]
clustercolors[is.na(clustercolors)] <- "gray80"
clustercolors[!is.sig] <- "gray80"
heatmapdata <- t(as.matrix(clustercentroids[,plottingChannels]))
labCols <- paste0("PhenoGraphMetacluster_",sprintf("%02d", clustercentroids$PhenoGraphMetacluster[]))
labCols[is.sig] <- paste0(labCols[is.sig],"*")
pdf(file = file.path(outdir2, paste0(analysisname," Rphenograph metacluster medians - plotting markers.pdf")), width=15, height=20)
heatmap.2(heatmapdata, col=mycolorfun, labCol=labCols, 
          labRow = key.channels$MarkerName[key.channels$MarkerChannel %in% plottingChannels], ColSideColors = clustercolors, colCol = clustercolors,
          trace="none", margins = c(25, 30), cexRow = 2, cexCol = 2, main=paste0(analysisname, "\n", "Rphenograph metacluster medians (plotting markers)"))
dev.off()


# Heatmap of metacluster medians, with significant metaclusters noted and rows scaled
pdf(file = file.path(outdir2, paste0(analysisname," Rphenograph metacluster medians - plotting markers - scaled rows.pdf")), width=15, height=20)
heatmap.2(heatmapdata, col=mycolorfun, scale="row", labCol=labCols, 
          labRow = key.channels$MarkerName[key.channels$MarkerChannel %in% plottingChannels], ColSideColors = clustercolors, colCol = clustercolors,
          trace="none", margins = c(25, 30), cexRow = 2, cexCol = 2, main=paste0(analysisname, "\n", "Rphenograph metacluster medians (plotting markers)"))
dev.off()

# Create boxplots of per-sample frequencies of all clusters
allclusters <- rownames(results)
plotdata <- cbind.data.frame(FCSprettyname=key.fcsids$FCSprettyname[has_condition_mask], conditions=conditions[has_condition_mask], t(as.data.frame(clusterpercentages)[allclusters,]))
rownames(plotdata) <- NULL
plotdata2 <- gather(data=plotdata, key=Metacluster, value=Percentage, allclusters)
plotdata2$Percentage[plotdata2$Percentage==0] <- 2^-10
plotdata2$Metacluster[plotdata2$Metacluster %in% sort(unique(plotdata2$Metacluster))[is.sig]] <- paste0(plotdata2$Metacluster[plotdata2$Metacluster %in% sort(unique(plotdata2$Metacluster))[is.sig]], "*")
pdf(file = file.path(outdir2, paste0(analysisname," Rphenograph boxplots of all metaclusters (raw values).pdf")), width=(length(allclusters)/5)+3, height=5)
p <- ggplot(plotdata2, aes(Metacluster, Percentage))
print(p + geom_boxplot(aes(fill = conditions, colour = conditions)) +
        geom_point(aes(Metacluster, Percentage, fill = conditions), position=position_dodge(width=0.75), size=0.5, show.legend=FALSE) +
        scale_colour_manual(values=c("grey30", "navy")) + scale_fill_manual(values=c("grey70", "#3977AF")) + 
        theme(panel.background = element_rect(fill="grey90"), axis.text.x = element_text(angle = 90, hjust = 1, colour = clustercolors)) + expand_limits(y=0) +
        ggtitle(analysisname, subtitle = paste("edgeR all metaclusters (raw scale)")))
dev.off()
pdf(file = file.path(outdir2, paste0(analysisname," Rphenograph boxplots of all metaclusters (log scale).pdf")), width=(length(allclusters)/5)+3, height=5)
p <- ggplot(plotdata2, aes(Metacluster, log2(Percentage)))
print(p + geom_boxplot(aes(fill = conditions, colour = conditions)) +
        geom_point(aes(Metacluster, log2(Percentage), fill = conditions), position=position_dodge(width=0.75), size=0.5, show.legend=FALSE) +
        scale_colour_manual(values=c("grey30", "navy")) + scale_fill_manual(values=c("grey70", "#3977AF")) + 
        theme(panel.background = element_rect(fill="grey90"), axis.text.x = element_text(angle = 90, hjust = 1, colour = clustercolors)) + expand_limits(y=0) +
        ggtitle(analysisname, subtitle = paste("edgeR all metaclusters (log scale)")))
dev.off()


# Create boxplots of per-sample frequencies of significant clusters
if(sum(is.sig) >= 1){
  sigclusters <- rownames(results)[is.sig]
  sigclustercolors <- clustercolors[as.numeric(sigclusters)]
  plotdata <- cbind.data.frame(FCSprettyname=key.fcsids$FCSprettyname[has_condition_mask], conditions=conditions[has_condition_mask], t(as.data.frame(clusterpercentages)[sigclusters,]))
  rownames(plotdata) <- NULL
  plotdata2 <- gather(data=plotdata, key=Metacluster, value=Percentage, sigclusters)
  plotdata2$Percentage[plotdata2$Percentage==0] <- 2^-10
  plotdata2$Metacluster <- paste0(plotdata2$Metacluster, "*")
  pdf(file = file.path(outdir2, paste0(analysisname," Rphenograph boxplots of significant metaclusters (raw values).pdf")), width=(length(sigclusters)/5)+3, height=5)
  p <- ggplot(plotdata2, aes(Metacluster, Percentage))
  print(p + geom_boxplot(aes(fill = conditions, colour = conditions)) +
          geom_point(aes(Metacluster, Percentage, fill = conditions), position=position_dodge(width=0.75), size=0.5, show.legend=FALSE) +
          scale_colour_manual(values=c("grey30", "navy")) + scale_fill_manual(values=c("grey70", "#3977AF")) +
          theme(panel.background = element_rect(fill="grey90"), axis.text.x = element_text(angle = 90, hjust = 1, colour = sigclustercolors)) + expand_limits(y=0) +
          ggtitle(analysisname, subtitle = paste("edgeR significant metaclusters (raw values)")))
  dev.off()
  pdf(file = file.path(outdir2, paste0(analysisname," Rphenograph boxplots of significant metaclusters (log scale).pdf")), width=(length(sigclusters)/5)+3, height=5)
  p <- ggplot(plotdata2, aes(Metacluster, log2(Percentage)))
  print(p + geom_boxplot(aes(fill = conditions, colour = conditions)) +
          geom_point(aes(Metacluster, log2(Percentage), fill = conditions), position=position_dodge(width=0.75), size=0.5, show.legend=FALSE) +
          scale_colour_manual(values=c("grey30", "navy")) + scale_fill_manual(values=c("grey70", "#3977AF")) + 
          theme(panel.background = element_rect(fill="grey90"), axis.text.x = element_text(angle = 90, hjust = 1, colour = sigclustercolors)) + expand_limits(y=0) +
          ggtitle(analysisname, subtitle = paste("edgeR significant metaclusters (log scale)")))
  dev.off()
} else {
  print(paste("[Warning] No clusters had uncorrected p-values below 0.05"))
}


# Save results
save.image(file = file.path(outdir2, paste0(analysisname," edgeR results.Rdata")))
