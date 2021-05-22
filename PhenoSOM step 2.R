# PhenoSOM of CD45-positive metaclusters
# Author:  Erin Simonds
# Version: 1.210521
#
# Summary: This script is meant to be run AFTER an initial round of PhenoSOM on the unclustered data. It
# takes PhenoGraph clusters as input, selects those that are CD45-positive, deconstructs them back
# to the original single-cell data, and then runs the PhenoSOM workflow on those cells. The reason
# for this sequential workflow is that the first round of PhenoSOM clustering does a good job of identifying
# putative immune cells (CD45+) and segregating them from tumor cells, but a poor job of identifying
#
# Some key parameters can be modified:
#
# 1. somx=10 and somy=10
# The default SOM grid is 10 x 10 = 100 nodes. This has worked well for me with datasets containing about
# 40 FCS files and about 2 million CD45-positive cells. I don't know how this would perform on much larger
# or smaller datasets. You can make the SOM map bigger or smaller by modifing xdim and ydim.
#
# 2. Rphenograph_k <- "k_30
# The default PhenoGraph K nearest neighbors constant "k" is 30. A lower value will result in more PhenoGraph clusters,
# which is not necessarily a good thing. I strive for as few PhenoGraph clusters as possible while still capturing
# the variety of cell types that exist in the dataset. I would recommend inspecting the PNG files produced by 
# this script that have filenames ending in "tSNE map of SOM nodes by Rphenograph metacluster k_xx.png" and choosing
# PhenoGraph "k" based on the value of "xx" in the filename. How do you pick the best one? This is subjective, but
# I look for the lowest "k" at which the SOM nodes fall into contiguous groups on the tSNE plot. There are other
# helpful files including "barplot of Rphenograph metacluster sizes (SOM nodes) with k_xx.png" and 
# "heatmap of sample contributions (SOM nodes) to Rphenograph metaclusters with k_xx.png" that can further inform
# this choice.
#
# 3. threshold.channel <- "Y89Di"
# You must set the channel (antibody) that will be used to identify metaclusters from the previous run of PhenoSOM
# that should be included in this script. If using the Fluidigm Y89-anti-CD45 antibody, the default setting of
# "Y89Di" is probably right.
#
# 4. threshold.value <- NULL
# This variable does not need to be set in advance of running the script, since it will ask the user for
# a threshold values upon the first run of the script. This value can be hard-coded to avoid the interactive
# dialog box in future runs. This variable decides the threshold for including metaclusters from the previous
# run of PhenoSOM. Only metaclusters with a median of the threshold.channel greater than this value will be
# included.
#
#
# Overall workflow of this script:
#     Import the old PhenoSOM metacluster results
#     Manually select a threshold value for threshold.channel; metaclusters with median intensity above this will be included
#     For each included metacluster, extract the single-cell data, keeping track of which FCS file each cell came from
#     Merge all of the single-cell data into a giant data.frame, with a column including the numerical GUID of the source file
#     Cluster each sample individually with FlowSOM (without metaclustering!), subsetting from the giant data.frame
#     Store the results of the FlowSOM runs in a list of FlowSOM objects
#     Compile the centroids of the FlowSOM nodes (sommap$codes) from each run into a single data.frame, including the numerical GUID of the source file
#     Run Rtsne on the data.frame of FlowSOM nodes
#     Save plot of tSNE map of SOM nodes by sample source
#     Run Rphenograph on the data.frame of FlowSOM nodes with the selected K (optional: You can run multiple Ks in a for-loop)
#     Plot Rphenograph metaclusters on Rtsne map with the selected K
#     Create a key linking FlowSOM nodes to Rphenograph metacluster assignments
#     Add the Rphenograph metacluster and FlowSOM node assignments as a column to the giant data.frame of cells
#     Aggregate the giant data.frame of single-cell data, calculating the mean of each marker in each Rphenograph metacluster
#     Aggregate the giant data.frame of single-cell data, calculating the median of each marker in each Rphenograph metacluster
#     Aggregate the giant data.frame of single-cell data, calculating the frequency (percentage) of each Rphenograph metacluster in each sample

rm(list = ls())

library(rstudioapi)
library(flowCore)
library(FlowSOM)
library(Rtsne)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(tidyr)
library(cytofkit2)
library(data.table)

# Set working directory to the location of the script, which has been Sourced from within Rstudio
setwd(dirname(rstudioapi::callFun("getActiveDocumentContext")$path)) # Set working directory to location of currently open Rstudio script
startdir <- getwd()
set.seed(42)

# Load previous clustering of all tumor cells (Rdata file) from Step 1
step1results <- dir("PhenoSOM_Step1_output", pattern="\\ results\\.Rdata")
load(step1results)

################################################# User-defined parameters ##################################################

# Select CD45 channel and cutoff
# Note:  This value can be hard-coded to avoid the dialog box in future runs)
threshold.channel <- "Y89Di"
threshold.value <- NULL  # Default: NULL This value can be hard-coded to avoid the interactive dialog box in future runs.

# Manually select an Rphenograph K to use for subsequent analysis
Rphenograph_k <- 30

# Define the size of the SOM map (see notes above for more details)
somx = 10
somy = 10

# Set a project name for this run of PhenoSOM on CD45pos clusters only
projectname <- "Your project name PhenoSOM Step 2"

outdir <- "PhenoSOM_Step2_output"
dir.create(outdir, showWarnings = F)

############################################################################################################################

# Log total time of script
print(paste("Script started at",Sys.time()))
timing <- as.data.frame(matrix(nrow=1, ncol=3, dimnames=list(c("script"), c("start", "stop", "elapsed"))))
timing["script","start"] <- Sys.time()
timing["script","stop"] <- Sys.time()
timing$elapsed <- timing$stop - timing$start

################################################ Declare functions #########################################################

# heatplot to colored by expression (scale range is 2nd to 98th percentile)
plotmap <- function(plotmapoutdir=outdir, marker, dotsize=0.5) {
  
  # Set color scale for plots
  mypal <- colorRampPalette(brewer.pal(9,"YlOrRd"))(256)
  
  # Load data for plot
  dat <- mapdata[,marker]
  
  # Label files and plots accordingly if the parameter was used for tSNE mapping
  if(marker %in% clustercolumns){
    titl <- paste0(prettynames[marker],"\nAll SOM nodes","\n(marker used to build tSNE map)")
  } else {
    titl <- paste0(prettynames[marker],"\nAll SOM nodes","\n(marker not used to build tSNE map")
  }
  pngname <- file.path(plotmapoutdir, paste0(projectname, " - ", prettynames[marker]," - All SOM nodes.png"))
  
  datmin <- quantile(dat,probs=0.02)
  datmax <- quantile(dat,probs=0.98)
  dat[dat<=datmin] <- datmin
  dat[dat>=datmax] <- datmax
  
  mycolors <- mypal[as.numeric(cut(dat,breaks = 256))]
  # Save plot to png
  png(file=pngname,width=1165,height=1024,units="px")
  
  par(mar=c(5.1, 4.1, 4.1, 2.1),bg="white",fg="black",col.axis="black",col.lab="black",col.main="black",col.sub="black")
  layout(matrix(1:2,ncol=2), width = c(5,1),height = c(1,1))
  plot(map, col=mycolors, pch=16, cex=dotsize, xlab="", ylab="", main=paste(titl,sep="\n"), sub=paste("min = ",round(datmin,2),"; max = ",round(datmax,2),sep=""))
  par(mar=c(6,1,6,3))
  image(1,seq(round(datmin,2),round(datmax,2),len=256),matrix(1:256,nrow=1),col=mypal,axes=FALSE,xlab="",ylab="")
  axis(2)
  
  dev.off()
  
  print(paste0(date()," - Printed ", basename(pngname)))
  
}

# heatplot to colored by expression (scale range is 2nd to 98th percentile) for specific samples
plotmap2 <- function(plotmapoutdir=outdir, marker, samplename, dotsize=0.5) {
  
  # Set color scale for plots
  mypal <- colorRampPalette(brewer.pal(9,"YlOrRd"))(256)
  
  # Load data for plot
  dat <- mapdata2[,marker]
  
  # Label files and plots accordingly if the parameter was used for tSNE mapping
  if(marker %in% clustercolumns){
    titl <- paste0(prettynames[marker],"\n",samplename,"\n(marker used to build tSNE map)")
  } else {
    titl <- paste0(prettynames[marker],"\n",samplename,"\n(marker not used to build tSNE map")
  }
  pngname <- file.path(plotmapoutdir, paste0(projectname, " - ", prettynames[marker], " - ", samplename, ".png"))
  
  datmin <- quantile(dat,probs=0.02)
  datmax <- quantile(dat,probs=0.98)
  dat[dat<=datmin] <- datmin
  dat[dat>=datmax] <- datmax
  
  mycolors <- mypal[as.numeric(cut(dat,breaks = 256))]
  # Save plot to png
  png(file=pngname,width=1165,height=1024,units="px")
  
  par(mar=c(5.1, 4.1, 4.1, 2.1),bg="white",fg="black",col.axis="black",col.lab="black",col.main="black",col.sub="black")
  layout(matrix(1:2,ncol=2), width = c(5,1),height = c(1,1))
  plot(map, col=mycolors, xlim=tsne_xlim, ylim=tsne_ylim, pch=16, cex=dotsize, xlab="", ylab="", main=paste(titl,sep="\n"), sub=paste("min = ",round(datmin,2),"; max = ",round(datmax,2),sep=""))
  par(mar=c(6,1,6,3))
  image(1,seq(round(datmin,2),round(datmax,2),len=256),matrix(1:256,nrow=1),col=mypal,axes=FALSE,xlab="",ylab="")
  axis(2)
  
  dev.off()
  
  print(paste0(date()," - Printed ", basename(pngname)))
  
}

# Read in all FCS files and apply asinh5 transform
timing["flowSet","start"] <- Sys.time()
print(paste("Reading in", length(dir(path=FCSdir, pattern=FCSfilenamepattern)), "FCS files started at",Sys.time()))
fcsorig <- read.flowSet(files=NULL, path=FCSdir, pattern=FCSfilenamepattern, alter.names = TRUE, transformation = FALSE, column.pattern = ".*[0-9]Di") # add which.lines=1:1000 for testing
timing["flowSet","stop"] <- Sys.time()
timing$elapsed <- timing$stop - timing$start
print(paste("Reading FCS files took", round(timing["flowSet",]$elapsed, digits=2), "seconds"))
numfcsfiles <- length(fcsorig)
print(paste("Loaded",numfcsfiles,"FCS files from",FCSdir))
asinhTrans <- arcsinhTransform(transformationId="cytofAsinh", a=0, b=0.2, c=0)
translist <- transformList(from=colnames(fcsorig[[1]]), tfun=asinhTrans)

print(paste("FCS file transformation started at",Sys.time()))
timing["fcstrans","start"] <- Sys.time()
fcsasinh <- transform(fcsorig, translist)
timing["fcstrans","stop"] <- Sys.time()
timing$elapsed <- timing$stop - timing$start
print(paste("Transforming FCS files took", round(timing["fcstrans",]$elapsed, digits=2), "seconds"))

# Define FCS numerical indices (FCSindex) and text GUIDs (FCSGUID)
FCSprettynames <- gsub(x=fsApply(fcsorig, FUN=function(x) identifier(x)), pattern="\\.fcs$", replacement="")
FCSsizes <- fsApply(fcsorig, nrow)
FCSindices <- as.matrix(1:length(fcsorig))
key.fcsids <- cbind.data.frame(FCSGUID=rownames(FCSprettynames), FCSprettyname=FCSprettynames, FCSindex=FCSindices, CellNum=FCSsizes)
write.csv(key.fcsids, file=file.path(outdir,paste0(projectname," sample IDs and sizes in input dataset.csv")), row.names = F)

# Collect a table of marker names in all FCS files
fcs.marker.names <- fsApply(fcsorig, FUN=function(x) pData(parameters(x))$desc)
write.csv(fcs.marker.names, file=file.path(outdir,paste0(projectname," original FCS marker names in input dataset.csv")), row.names = T)

# Create key table linking FCS channels to names using last FCS file in dataset
key.channels <- cbind.data.frame(MarkerIndex=(1:nrow(pData(parameters(fcsorig[[length(fcsorig)]])))), MarkerChannel=pData(parameters(fcsorig[[length(fcsorig)]]))$name, MarkerName=pData(parameters(fcsorig[[length(fcsorig)]]))$desc)
key.channels$MarkerPrettyName <- paste0(key.channels$MarkerChannel,"<",key.channels$MarkerName,">")
key.channels$IsClusteringChannel <- as.numeric(key.channels$MarkerChannel %in% clusteringChannels)
key.channels$IsPlottingChannel <- as.numeric(key.channels$MarkerChannel %in% plottingChannels)
row.names(key.channels) <- key.channels$MarkerChannel
write.csv(key.channels, file=file.path(outdir,paste0(projectname," channel and marker names used.csv")), row.names = F)

# Merge all single-cell data into giant data.frame
print(paste("Merging all cells into giant data frame started at",Sys.time()))
timing["merging","start"] <- Sys.time()
cells.list <- fsApply(fcsasinh, FUN=function(x) cbind.data.frame(FCSindex=key.fcsids[identifier(x),"FCSindex"], CellIndex=1:nrow(exprs(x)), exprs(x)))
cells.df <- do.call(rbind.data.frame, cells.list)
rownames(cells.df) <- NULL
timing["merging","stop"] <- Sys.time()
timing$elapsed <- timing$stop - timing$start
print(paste("Merging all cells into giant data frame took", round(timing["merging",]$elapsed/60, digits=2), "minutes"))

# Add the Rphenograph metacluster and FlowSOM node assignments as a column to the giant data.frame of cells
som.mapping.list <- lapply(X=sommaps.list, FUN=function(x) {cbind.data.frame(FCSindex=x$FCSindex, CellIndex=1:nrow(x$mapping), SOMnode=x$mapping[,1])})
som.mapping.df <- do.call(rbind.data.frame, som.mapping.list)
clustermapping.df <- merge(som.mapping.df, key.clusters)
clustermapping.df <- clustermapping.df[order(clustermapping.df$FCSindex, clustermapping.df$CellIndex, clustermapping.df$SOMnode),]
cellswithclusters.df <- merge(clustermapping.df, cells.df)


############# Select cutoff #################

# Select CD45 cutoff
threshold.channel <- "Y89Di"
threshold.value <- 0.4  # this value can be hard-coded to avoid the interactive dialog box in future runs

# If threshold.value is undefined, make a dot plot of original metacluster medians and ask the user to choose the gating threshold
if(is.null(threshold.value) | length(threshold.value)==0){
  p <- qplot(x=sprintf("%02d", phenograph.metaclusters.exprs.median.df$PhenoGraphMetacluster), y=phenograph.metaclusters.exprs.median.df[,threshold.channel], xlab="PhenoSOM Metacluster", ylab=paste0("Median of ",threshold.channel), ymin=0)
  thresholdplotfilename <- file.path(outdir, paste0(projectname," plot of original PhenoSOM metacluster medians.png"))
  ggsave(thresholdplotfilename, p, width = 8, height = 8)  
}

# Save a dot plot of original metacluster medians with the gating threshold overlaid
outdir7 <- file.path(outdir, "Plots of original PhenoSOM clusters")
dir.create(outdir7, showWarnings = F)
png(filename = file.path(outdir7, paste0(projectname," plot of original PhenoSOM metacluster medians with gating threshold.png")), width = 1000, height = 200)
  p <- ggplot(data=phenograph.metaclusters.exprs.median.df, aes_string(x="PhenoGraphMetacluster", y=threshold.channel)) + geom_point() + xlab("PhenoSOM Metacluster") + ylab(paste0("Median of ",threshold.channel))
  print(p)
dev.off()

if(is.null(threshold.value) | length(threshold.value)==0){
  rstudioapi::showDialog(title="Choose CD45 threshold", message=paste0("View the file '",thresholdplotfilename,"' and choose a threshold based on the plot. Metaclusters with median intensities above this value will be included. Define the threshold in the next prompt."))
  threshold.value <- as.numeric(rstudioapi::showPrompt(title="Choose CD45 threshold", message=paste0("Inspect the file ",thresholdplotfilename," and select a threshold based on the plot. Metaclusters with median intensities above this value will be included."), default="1.0"))
}

if(is.null(threshold.value) | length(threshold.value)==0){
  stop("You must set threshold.value to a number above 0 to continue")    
}

# Keep metaclusters above the threshold value
metaclusters.keep.mask <- phenograph.metaclusters.exprs.median.df[,threshold.channel] > threshold.value
metaclusters.keep <- sort(as.integer(phenograph.metaclusters.exprs.median.df$PhenoGraphMetacluster)[metaclusters.keep.mask])

# Save a dot plot of original metacluster medians with the gating threshold overlaid
png(filename = file.path(outdir7, paste0(projectname," plot of original PhenoSOM metacluster medians with ", threshold.channel," gating threshold.png")), width = 1000, height = 200)
p <- ggplot(data=phenograph.metaclusters.exprs.median.df, aes_string(x="PhenoGraphMetacluster", y=threshold.channel, color=metaclusters.keep.mask)) + geom_point() + xlab("PhenoSOM Metacluster") + ylab(paste0("Median of ",threshold.channel)) + 
  geom_hline(yintercept = threshold.value, colour="red") + scale_color_discrete(name="Keep")
print(p)
dev.off()

# Save a violin plot of original metaclusters with the gating threshold overlaid
plotdata <- cbind.data.frame(PhenoGraphMetacluster=sprintf("%02d", cellswithclusters.df$PhenoGraphMetacluster), cellswithclusters.df[,threshold.channel])
colnames(plotdata)[2] <- threshold.channel
plotdata <- cbind(plotdata, Keep=plotdata$PhenoGraphMetacluster %in% sprintf("%02d", metaclusters.keep))
png(filename = file.path(outdir7, paste0(projectname," plot of original PhenoSOM metacluster violins with ", threshold.channel," gating threshold.png")), width = 1000, height = 200)
p <- ggplot(plotdata, aes_string(x="PhenoGraphMetacluster", y=threshold.channel, fill="Keep"))
print(p + geom_violin(aes(group=PhenoGraphMetacluster), scale="width") + geom_hline(yintercept = threshold.value, colour="red"))
dev.off()

############# Restrict the rest of the analysis to the metaclusters that are above the cutoff #################
cellswithclusters.df.restricted <- cellswithclusters.df[cellswithclusters.df$PhenoGraphMetacluster %in% metaclusters.keep, ]

# Fix CellIndex values in cells.df.restricted so they are incrementing consistently (this becomes important later in SOM mapping step)
cells.df <- cellswithclusters.df.restricted[ , !colnames(cellswithclusters.df) %in% c("SOMnode", "PhenoGraphMetacluster")]
cells.df.2 <- cells.df[order(cells.df$FCSindex, cells.df$CellIndex),]
rownames(cells.df.2) <- NULL
for(fi in unique(cells.df.2$FCSindex)){
  cells.df.2[cells.df.2$FCSindex==fi,"CellIndex"] <- 1:nrow(cells.df.2[cells.df.2$FCSindex==fi,])
}
cells.df <- cells.df.2

suppressWarnings(rm(allobjects))
suppressWarnings(rm(removeobjects))
allobjects <- ls(all.names=TRUE)
removeobjects <- allobjects[!((allobjects %in% "cells.df") | 
                                (allobjects %in% "plotmap") | 
                                (allobjects %in% "plotmap2") | 
                                (allobjects %in% "timing") | 
                                (allobjects %in% "startdir") | 
                                (allobjects %in% "key.fcsids") | 
                                (allobjects %in% "key.channels") |  
                                (allobjects %in% "FCSdir") |   
                                (allobjects %in% "FCSindices") |
                                (allobjects %in% "FCSfilenamepattern") | 
                                (allobjects %in% "projectname") | 
                                (allobjects %in% "outdir") | 
                                (allobjects %in% "FCSprettynames") |
                                (allobjects %in% "threshold.channel") |
                                (allobjects %in% "threshold.value") |
                                (allobjects %in% "Rphenograph_k") |
                                (allobjects %in% "somx") |
                                (allobjects %in% "somy") |
                                (allobjects %in% "clusteringMarkers") |
                                (allobjects %in% "clusteringChannels") |
                                (allobjects %in% "plottingMarkers") |
                                (allobjects %in% "plottingChannels") |
                                (allobjects %in% ".Random.seed"))]
rm(list = removeobjects)

# Set markers to use for clustering
if(length(dir(getwd(), pattern=".*\\ clustering markers\\.txt"))!=1) print("[Error] clustering markers.txt file not found.")
clusteringMarkers <- as.character(unlist(read.delim(file=dir(getwd(), pattern=".*\\ clustering markers\\.txt"), header=F)))
clusteringChannels <- gsub(x=clusteringMarkers, pattern="<.*>", replacement="")
print(paste("Markers used for clustering:"))
print(as.matrix(clusteringMarkers))

# Set markers to use for plotting
if(length(dir(getwd(), pattern=".*\\ plotting markers\\.txt"))!=1) print("[Error] plotting markers.txt file not found.")
plottingMarkers <- as.character(unlist(read.delim(file=dir(getwd(), pattern=".*\\ plotting markers\\.txt"), header=F)))
plottingChannels <- gsub(x=plottingMarkers, pattern="<.*>", replacement="")
print(paste("Markers used for plotting:"))
print(as.matrix(plottingMarkers))

# Define FCS numerical indices (FCSindex) and text GUIDs (FCSGUID)
key.fcsids$CellNum <- table(cells.df$FCSindex)
write.csv(key.fcsids, file=file.path(outdir, paste0(projectname," sample IDs and sizes in input dataset.csv")), row.names = F)

# Check if Rdata file from previous run of PhenoSOM exists, and load it if present. Otherwise, run fresh.
if(!exists("outdir") || outdir==""){
  print("Folder not specified. Try again.")
  stop()
} else {
  if(sum(grepl(pattern=".*PhenoSOM.*results\\.Rdata$", x=dir(".", full.names=F)))<1){
    print("Previous PhenoSOM result not found. Running clustering from scratch.")
    newrun <- TRUE
  } else {
    if(sum(grepl(pattern=".*PhenoSOM.*results\\.Rdata$", x=dir(".", full.names=F)))>1){
      print("More than one previous PhenoSOM result found. Try again.")
      stop()
    } else {
      print("Previous PhenoSOM result found.")
      print(paste0("Loading previous Rdata file:  ", grep(pattern=".*PhenoSOM.*results\\.Rdata$", x=dir(".", full.names=F), value=T)))
      load(grep(pattern=".*PhenoSOM.*results\\.Rdata$", x=dir(".", full.names=F), value=T))
      newrun <- FALSE
    }
  }
}


# Run clustering if Rdata file was not loaded
if(newrun){
  # Run FlowSOM on each file separately and store as a list
  timing["totalclust","start"] <- Sys.time()
  print(paste("SOM node assignment (no metaclustering) started at ", timing["totalclust","start"] ))
  sommaps.list <- list()
  for(f1 in FCSindices){
    
    # Load expression data
    print("----------------------------------------------")
    print(paste0("Processing file:  ", key.fcsids$FCSGUID[f1]))
    clustdata <- as.matrix(subset(cells.df, FCSindex==f1)[,clusteringChannels])
    print(paste("Clustering on", nrow(clustdata),"cells using",ncol(clustdata),"out of",ncol(cells.df),"total channels"))
    colnames(clustdata) <- key.channels[clusteringChannels,"MarkerName"]
    
    # Build initial clustering SOM matrix
    # This takes 19 seconds for 60000 cells and 31 clustering markers with k=50
    print(paste("Building SOM started at", date()))
    som.start <- Sys.time()
    sommap <- SOM(clustdata, silent = F, xdim=somx, ydim=somy)
    sommap$FCSindex <- f1
    print(paste("Building SOM ended at", date()))
    som.end <- Sys.time()
    som.elapsed <- (som.end - som.start)
    print(paste("Building SOM took", round(as.numeric(som.elapsed, units="secs"), digits=2), "seconds"))
    
    # Store the results of all FlowSOM runs in a list of FlowSOM objects
    sommaps.list[[as.character(key.fcsids$FCSGUID[f1])]] <- sommap
    
  }
  timing["totalclust","stop"] <- Sys.time()
  timing$elapsed <- timing$stop - timing$start
  print("**********************************************")
  print(paste("FlowSOM clustering all files took", round(timing["totalclust",]$elapsed/60, digits=2), "minutes"))
  print("**********************************************")
  
  # Compile the centroids of the FlowSOM nodes (sommap$codes) from each run into a single data.frame, including the numerical GUID of the source file
  som.centroids.list <- lapply(X=sommaps.list, FUN=function(x) {cbind.data.frame(FCSindex=x$FCSindex, SOMnode=1:nrow(x$codes), x$codes)})
  som.centroids.df <- do.call(rbind.data.frame, som.centroids.list)
  
  # Run Rtsne on the data.frame of FlowSOM nodes
  print(paste("Rtsne of FlowSOM nodes started at",Sys.time()))
  timing["tsne","start"] <- Sys.time()
  tsnemap <- Rtsne(som.centroids.df[,-c(1:2)], verbose=T)
  timing["tsne","stop"] <- Sys.time()
  timing$elapsed <- timing$stop - timing$start
  print(paste("Rtsne took", round(timing["tsne",]$elapsed, digits=2), "seconds"))
  

  # Run Rphenograph on the data.frame of FlowSOM nodes with various Ks and show metaclusters on Rtsne map
  phenograph.metaclusters.list <- list()
  print(paste("Rphenograph of FlowSOM nodes started at",Sys.time()))
  
  Rphenograph_k_ID <- paste0("k_", Rphenograph_k)
  for(k in c(Rphenograph_k)){  #OPTIONAL:  A range of values of K can be set as a vector here
  	timing["phenograph","start"] <- Sys.time()
  	phenograph.start <- Sys.time()
  	phenograph.metaclusters.list[[Rphenograph_k_ID]] <- as.numeric(membership(Rphenograph(som.centroids.df[,-c(1:2)], k=Rphenograph_k)))
  	timing["phenograph","stop"] <- Sys.time()
  	timing$elapsed <- timing$stop - timing$start
  	print(paste("R phenograph with k =", k, "took", round(timing["phenograph",]$elapsed, digits=2), "seconds"))
  }
}

# Save plot of tSNE map of SOM nodes by sample source
outdir4 <- file.path(outdir, "Maps and plots of SOM nodes")
dir.create(outdir4, showWarnings = F)
png(filename = file.path(outdir4,paste0(projectname," tSNE map of SOM nodes by sample source.png")), width = 1000, height = 1000)
plotdata <- cbind.data.frame(tsnemap$Y, som.centroids.df$FCSindex)
colnames(plotdata) <- c("tSNE_1", "tSNE_2", "FCSindex")
mycolorfun <- colorRampPalette(brewer.pal(11,"Spectral"))
mycolors <- rev(mycolorfun(nrow(key.fcsids)))
p <- ggplot(plotdata, aes(tSNE_1, tSNE_2))
print(p + geom_point(aes(colour = as.factor(FCSindex))) +
        scale_colour_manual("FCSindex", values=mycolors, guide="legend") +
        ggtitle(paste(projectname, "SOM nodes by sample source"), subtitle = NULL))
dev.off()

# Plot SOM nodes with Rphenograph metaclusters overlaid various Ks
for(Rpk in names(phenograph.metaclusters.list)){
  print(paste("Plotting phenograph metaclusters for",Rpk))
  png(filename = file.path(outdir4,paste0( projectname, " tSNE map of SOM nodes by Rphenograph metacluster ", Rpk, ".png")), width = 1000, height = 1000)
  plotdata <- cbind.data.frame(tsnemap$Y, som.centroids.df$FCSindex, phenograph.metaclusters.list[[Rpk]])
  colnames(plotdata) <- c("tSNE_1", "tSNE_2", "FCSindex", "PhenoGraphMetacluster")
  mycolors <- rainbow(length(unique(plotdata$PhenoGraphMetacluster)))
  p <- ggplot(plotdata, aes(tSNE_1, tSNE_2))
  print(p + geom_point(aes(colour = as.factor(PhenoGraphMetacluster))) +
          scale_colour_manual("Cluster ID", values=mycolors, guide="legend") +
          ggtitle(paste(projectname, "SOM nodes by Rphenograph metacluster with",Rpk), subtitle = NULL))
  dev.off()
}

# Save barplots of Rphenograph metacluster sizes (SOM nodes)
mycolorfun <- colorRampPalette(brewer.pal(11,"Spectral"))
mycolors <- rev(mycolorfun(nrow(key.fcsids)))
for(Rpk in names(phenograph.metaclusters.list)){
  kc <- cbind.data.frame(FCSindex=som.centroids.df$FCSindex, SOMnode=som.centroids.df$SOMnode, PhenoGraphMetacluster=phenograph.metaclusters.list[[Rpk]])
  png(filename = file.path(outdir4,paste0(projectname," barplot of Rphenograph metacluster sizes (SOM nodes) with ",Rpk,".png")), width = 1000, height = 200)
  p <- ggplot(kc, aes(x = sprintf("%02d", PhenoGraphMetacluster), fill = sprintf("%02d", FCSindex)))
  print(p + geom_bar() + xlab("PhenoSOM Metacluster") + ylab("# of SOM nodes") + scale_fill_manual("FCSindex", values=mycolors, guide="legend") + guides(fill=guide_legend(nrow=1, label.position = "bottom")) + theme(legend.position="bottom") + ggtitle(paste0("Rphenograph metacluster sizes (SOM nodes) with ",Rpk)))
  dev.off()
}

# Save heatmaps of Rphenograph metacluster sizes (SOM nodes) colored by sample contribution
for(Rpk in names(phenograph.metaclusters.list)){
  kc <- cbind.data.frame(FCSindex=som.centroids.df$FCSindex, SOMnode=som.centroids.df$SOMnode, PhenoGraphMetacluster=phenograph.metaclusters.list[[Rpk]])
  png(filename = file.path(outdir4,paste0(projectname," heatmap of sample contributions (SOM nodes) to Rphenograph metaclusters with ",Rpk,".png")), width = 800, height = 500)
  SOMnodes.per.phenographmetacluster <- aggregate(SOMnode ~ PhenoGraphMetacluster + FCSindex, data=kc, FUN=length)
  colnames(SOMnodes.per.phenographmetacluster)[3] <- "NumSOMnodes"
  p <- ggplot(SOMnodes.per.phenographmetacluster)
  print(p + geom_tile(aes(x=sprintf("%02d", PhenoGraphMetacluster), y=sprintf("%02d", FCSindex), fill=NumSOMnodes)) + scale_fill_distiller("# of SOM nodes", palette = "Spectral") + xlab("PhenoSOM Metacluster") + ylab("FCS index"))
  dev.off()
}

# Create a key linking FlowSOM nodes to Rphenograph metacluster assignments
key.clusters <- cbind.data.frame(FCSindex=som.centroids.df$FCSindex, SOMnode=som.centroids.df$SOMnode, PhenoGraphMetacluster=phenograph.metaclusters.list[[Rphenograph_k_ID]])

# Add the Rphenograph metacluster and FlowSOM node assignments as a column to the giant data.frame of cells
som.mapping.list <- lapply(X=sommaps.list, FUN=function(x) {cbind.data.frame(FCSindex=x$FCSindex, CellIndex=1:nrow(x$mapping), SOMnode=x$mapping[,1])})
som.mapping.df <- do.call(rbind.data.frame, som.mapping.list)
clustermapping.df <- merge(som.mapping.df, key.clusters)
clustermapping.df <- clustermapping.df[order(clustermapping.df$FCSindex, clustermapping.df$CellIndex, clustermapping.df$SOMnode),]
cellswithclusters.df <- merge(clustermapping.df, cells.df)

# Save barplot of Rphenograph metacluster sizes (cell number) colored by sample contribution
outdir5 <- file.path(outdir, "Maps and plots of PhenoGraph clusters")
dir.create(outdir5, showWarnings = F)
mycolorfun <- colorRampPalette(brewer.pal(11,"Spectral"))
mycolors <- rev(mycolorfun(nrow(key.fcsids)))
png(filename = file.path(outdir5,paste0(projectname," barplot of Rphenograph metacluster sizes (cell number) with ",Rphenograph_k_ID,".png")), width = 1000, height = 200)
p <- ggplot(cellswithclusters.df, aes(x = sprintf("%02d", PhenoGraphMetacluster), fill = sprintf("%02d", FCSindex)))
print(p + geom_bar() + xlab("PhenoSOM Metacluster") + ylab("Cell count") + scale_fill_manual("FCSindex", values=mycolors, guide="legend") + guides(fill=guide_legend(nrow=1, label.position = "bottom")) + theme(legend.position="bottom") + ggtitle(paste0("Rphenograph metacluster sizes (cell number) with ",Rphenograph_k_ID)))
dev.off()

# Aggregating the giant data.frame of single-cell data, calculate the mean of each marker in each Rphenograph metacluster
print(paste("Calculating the mean of each marker in each Rphenograph metacluster started at",Sys.time()))
cellswithclusters.dt <- as.data.table(cellswithclusters.df)
selectedCols <- colnames(cellswithclusters.df)[-c(1:4)] # select columns to calculate means and medians
phenograph.metaclusters.exprs.mean.dt <- cellswithclusters.dt[, lapply(.SD, mean, na.rm=TRUE), by=PhenoGraphMetacluster, .SDcols=selectedCols]
setorder(phenograph.metaclusters.exprs.mean.dt, PhenoGraphMetacluster)
phenograph.metaclusters.exprs.mean.df <- as.data.frame(phenograph.metaclusters.exprs.mean.dt)

# Aggregating the giant data.frame of single-cell data, calculate the median of each marker in each Rphenograph metacluster
print(paste("Calculating the median of each marker in each Rphenograph metacluster started at",Sys.time()))
phenograph.metaclusters.exprs.median.dt <- cellswithclusters.dt[, lapply(.SD, median, na.rm=TRUE), by=PhenoGraphMetacluster, .SDcols=selectedCols]
setorder(phenograph.metaclusters.exprs.median.dt, PhenoGraphMetacluster)
phenograph.metaclusters.exprs.median.df <- as.data.frame(phenograph.metaclusters.exprs.median.dt)

# Subsetting the giant data.frame of single-cell data, calculate the frequency (percentage) of each SOM node in each sample
SOMnodes.counts.long.df <- aggregate(CellIndex ~ SOMnode + FCSindex, data=cellswithclusters.df[c("FCSindex","CellIndex","SOMnode")], FUN=length)
colnames(SOMnodes.counts.long.df)[3] <- "NumCells"
SOMnodes.counts.df <- SOMnodes.counts.long.df %>% spread(key=c(SOMnode), value=NumCells, fill=0)
colnames(SOMnodes.counts.df) <- c("FCSindex", paste0("SOMnode_",colnames(SOMnodes.counts.df[,-1])))
write.csv(SOMnodes.counts.df, file=file.path(outdir4,paste0(projectname," SOM node sizes (cell number) per file.csv")), row.names = F)

# calculate SOM node percentages on a per-sample basis
allFCScelltotals <- key.fcsids$CellNum
SOMnodes.percentages.df <- cbind.data.frame(FCSindex=key.fcsids$FCSindex, apply(as.matrix(SOMnodes.counts.df[,-1]), 2, function(x) (x/allFCScelltotals)*100))
write.csv(SOMnodes.percentages.df, file=file.path(outdir,paste0(projectname," SOM node sizes (percentage) per file.csv")), row.names = F)

# Subsetting the giant data.frame of single-cell data, calculate the frequency (percentage) of each Rphenograph metacluster in each sample
phenograph.metaclusters.counts.long.df <- aggregate(CellIndex ~ PhenoGraphMetacluster + FCSindex, data=cellswithclusters.df[c("FCSindex","CellIndex","PhenoGraphMetacluster")], FUN=length)
colnames(phenograph.metaclusters.counts.long.df)[3] <- "NumCells"
phenograph.metaclusters.counts.df <- phenograph.metaclusters.counts.long.df %>% spread(key=c(PhenoGraphMetacluster), value=NumCells, fill=0)
colnames(phenograph.metaclusters.counts.df) <- c("FCSindex", paste0("PhenoGraphMetacluster_",colnames(phenograph.metaclusters.counts.df[,-1])))
write.csv(phenograph.metaclusters.counts.df, file=file.path(outdir5,paste0(projectname," Rphenograph metacluster sizes (cell number) per file.csv")), row.names = F)

# calculate phenograph metacluster percentages on a per-sample basis
allFCScelltotals <- key.fcsids$CellNum
phenograph.metaclusters.percentages.df <- cbind.data.frame(FCSindex=key.fcsids$FCSindex, apply(as.matrix(phenograph.metaclusters.counts.df[,-1]), 2, function(x) (x/allFCScelltotals)*100))
write.csv(phenograph.metaclusters.percentages.df, file=file.path(outdir5,paste0(projectname," Rphenograph metacluster sizes (percentage) per file.csv")), row.names = F)

# Plot SOM nodes with Rphenograph metaclusters overlaid and scaled by size (cell number)
print(paste("Plotting phenograph metaclusters for",Rphenograph_k_ID))
png(filename = file.path(outdir4,paste0( projectname, " tSNE map of SOM nodes by Rphenograph metacluster ", Rphenograph_k_ID, " and size scaled to cell number.png")), width = 1000, height = 1000)
plotdata <- cbind.data.frame(tsnemap$Y, som.centroids.df$FCSindex, som.centroids.df$SOMnode, phenograph.metaclusters.list[[Rphenograph_k_ID]])
colnames(plotdata) <- c("tSNE_1", "tSNE_2", "FCSindex", "SOMnode", "PhenoGraphMetacluster")
plotdata2 <- merge(plotdata, SOMnodes.counts.long.df)
mycolors <- rainbow(length(unique(plotdata$PhenoGraphMetacluster)))
p <- ggplot(plotdata2, aes(tSNE_1, tSNE_2))
print(p + geom_point(aes(colour = as.factor(PhenoGraphMetacluster), size = NumCells)) +
        scale_colour_manual("Cluster ID", values=mycolors, guide="legend") +
        scale_size_area("Number of cells", max_size = 8) +
        ggtitle(paste(projectname, "SOM nodes by Rphenograph metacluster with",Rphenograph_k_ID,"and size scaled by cell number"), subtitle = NULL))
dev.off()

# make heatmap of phenograph metacluster percentages on a per-sample basis
mycolorfun <- colorRampPalette(rev(brewer.pal(11,"RdBu")))
pdf(file = file.path(outdir5,paste0(projectname," Rphenograph metacluster sizes (percentage) per file.pdf")), width=20, height=20)
heatmap.2(as.matrix(phenograph.metaclusters.percentages.df[,-1]), col=mycolorfun, labRow = key.fcsids$FCSprettyname, trace="none", margins = c(15, 40), cexCol = 1.2, main="Rphenograph metacluster sizes (percentage) per file")
dev.off()

# make heatmap of phenograph metacluster marker medians - plotting markers
mycolorfun <- colorRampPalette(rev(brewer.pal(11,"RdBu")))
pdf(file = file.path(outdir5,paste0(projectname," Rphenograph metacluster medians - plotting markers.pdf")), width=15, height=20)
heatmap.2(t(as.matrix(phenograph.metaclusters.exprs.median.df[,plottingChannels])), col=mycolorfun, labCol=paste0("PhenoGraphMetacluster_",phenograph.metaclusters.exprs.median.df$PhenoGraphMetacluster), labRow = key.channels$MarkerName[key.channels$MarkerChannel %in% plottingChannels], trace="none", margins = c(15, 20), cexRow = 1.2, cexCol = 1.2, main="Rphenograph metacluster medians (plotting markers)")
dev.off()

# make heatmap of phenograph metacluster marker medians - clustering markers only
mycolorfun <- colorRampPalette(rev(brewer.pal(11,"RdBu")))
pdf(file = file.path(outdir5,paste0(projectname," Rphenograph metacluster medians - clustering markers.pdf")), width=15, height=15)
heatmap.2(t(as.matrix(phenograph.metaclusters.exprs.median.df[,clusteringChannels])), col=mycolorfun, labCol=paste0("PhenoGraphMetacluster_",phenograph.metaclusters.exprs.median.df$PhenoGraphMetacluster), labRow = key.channels$MarkerName[key.channels$MarkerChannel %in% clusteringChannels], trace="none", margins = c(15, 20), cexRow = 1.2, cexCol = 1.2, main="Rphenograph metacluster medians (clustering markers only)")
dev.off()

# make heatmap of phenograph metacluster marker medians - plotting markers scaled to row
mycolorfun <- colorRampPalette(rev(brewer.pal(11,"RdBu")))
pdf(file = file.path(outdir5,paste0(projectname," Rphenograph metacluster medians - plotting markers (scaled).pdf")), width=15, height=20)
heatmap.2(t(as.matrix(phenograph.metaclusters.exprs.median.df[,plottingChannels])), scale=("row"), col=mycolorfun, labCol=paste0("PhenoGraphMetacluster_",phenograph.metaclusters.exprs.median.df$PhenoGraphMetacluster), labRow = key.channels$MarkerName[key.channels$MarkerChannel %in% plottingChannels], trace="none", margins = c(15, 20), cexRow = 1.2, cexCol = 1.2, main="Rphenograph metacluster medians (plotting markers)")
dev.off()

# make heatmap of phenograph metacluster marker medians - clustering markers only scaled to row
mycolorfun <- colorRampPalette(rev(brewer.pal(11,"RdBu")))
pdf(file = file.path(outdir5,paste0(projectname," Rphenograph metacluster medians - clustering markers (scaled).pdf")), width=15, height=15)
heatmap.2(t(as.matrix(phenograph.metaclusters.exprs.median.df[,clusteringChannels])), scale=("row"), col=mycolorfun, labCol=paste0("PhenoGraphMetacluster_",phenograph.metaclusters.exprs.median.df$PhenoGraphMetacluster), labRow = key.channels$MarkerName[key.channels$MarkerChannel %in% clusteringChannels], trace="none", margins = c(15, 20), cexRow = 1.2, cexCol = 1.2, main="Rphenograph metacluster medians (clustering markers only)")
dev.off()

# save CSV of Phenograph metacluster marker medians
phenograph.metaclusters.exprs.median.df.out <- phenograph.metaclusters.exprs.median.df
colnames(phenograph.metaclusters.exprs.median.df.out)[-1] <- key.channels$MarkerPrettyName
write.csv(phenograph.metaclusters.exprs.median.df.out, file=file.path(outdir5,paste0(projectname," Rphenograph metacluster medians - all markers.csv")), row.names = F)

# Save tSNE plots of SOM nodes colored by marker expression
mapdata <- som.centroids.df[,-c(1:2)]
clustermarkers <- as.character(subset(key.channels, IsClusteringChannel==1)$MarkerName)
clustercolumns <- 1:ncol(mapdata)
channelnames <- as.character(subset(key.channels, IsClusteringChannel==1)$MarkerChannel)
channeldescs <- gsub(x=gsub(x=gsub(clustermarkers, pattern=".*<", replacement=""), pattern=">", replacement=""), pattern="_", replacement=" ")
prettynames <- paste0(channeldescs, " - ", channelnames)
map <- tsnemap$Y
outdir2 <- file.path(outdir,"tSNE of SOM nodes - All SOM nodes by marker")
dir.create(outdir2, showWarnings = FALSE)
lapply(clustercolumns, plotmap, plotmapoutdir=outdir2, dotsize=1)

# Save per-sample tSNE plots of SOM nodes colored by marker expression
outdir3 <- file.path(outdir,"tSNE of SOM nodes - Each sample by marker")
dir.create(outdir3, showWarnings = FALSE)
for(ind in FCSindices){
  mapdata2 <- subset(som.centroids.df, FCSindex==ind)[,-c(1:2)]
  samplename <- FCSprettynames[ind]
  clustermarkers <- as.character(subset(key.channels, IsClusteringChannel==1)$MarkerName)
  clustercolumns <- 1:ncol(mapdata)
  channelnames <- as.character(subset(key.channels, IsClusteringChannel==1)$MarkerChannel)
  channeldescs <- gsub(x=gsub(x=gsub(clustermarkers, pattern=".*<", replacement=""), pattern=">", replacement=""), pattern="_", replacement=" ")
  prettynames <- paste0(channeldescs, " - ", channelnames)
  map <- tsnemap$Y[som.centroids.df$FCSindex %in% ind, ]
  bigmap <- tsnemap$Y
  tsne_xlim <- c(min(bigmap[,1]), max(bigmap[,1]))
  tsne_ylim <- c(min(bigmap[,2]), max(bigmap[,2]))
  outdir6 <- file.path(outdir3, samplename)
  dir.create(outdir6, showWarnings = F)
  Sys.sleep(1)
  lapply(clustercolumns, FUN=plotmap2, samplename=samplename, plotmapoutdir=outdir6, dotsize=2)
}

# Save R workspace to file except cells.df and cellswithclusters.df which are too huge
suppressWarnings(rm(allobjects))
suppressWarnings(rm(mostobjects))
allobjects <- ls(all.names=TRUE)
mostobjects <- allobjects[!((allobjects %in% "cells.df") | 
                              (allobjects %in% "cells.list") | 
                              (allobjects %in% "cellswithclusters.df") | 
                              (allobjects %in% "fcsorig") | 
                              (allobjects %in% "fcsasinh"))]
save(file=file.path(outdir,paste0(projectname," results.Rdata")), list=mostobjects, envir=.GlobalEnv)
print(paste("Workspace image saved at",Sys.time()))

# Log total time of script
print(paste("Script finished at",Sys.time()))

timing["script","stop"] <- Sys.time()
timing$elapsed <- timing$stop - timing$start
print("**********************************************")
print(paste("Script total run time:", round(timing["script",]$elapsed/60, digits=2), "minutes"))
print("**********************************************")
