# ==============================================================================
#      Annotate G4 and PRDM9 binding sites to genomic features
# ==============================================================================
# Author: Dr. Angelika Lahnsteiner
# Date: 2026-07-30
# ==============================================================================


# ------------------------------
# Load libraries
# ------------------------------
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(ggplot2)

# ------------------------------
# Set directory and load files
# ------------------------------
setwd("~/Desktop/G4_PRDM9_Analysis/Analysis/ChIPseeker")
samplefiles <- list.files("~/Desktop/G4_PRDM9_Analysis/Analysis/ChIPseeker/files", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
print(samplefiles)
names(samplefiles) <- c("G4s","PRDM9")
print(samplefiles)


# ------------------------------
# Annotate peaks to genomic features and prepare plots
# ------------------------------
# perform the next step for G4s and PRDM9 peaks seperately
peak <- readPeakFile(samplefiles[["G4s"]])
peak <- readPeakFile(samplefiles[["PRDM9"]])
peak

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

# plot peaks
plotAnnoPie(peakAnno) # pie plot
plotAnnoBar(peakAnnoList) # bar plot

# plot annotated G4 and PRDM9 peaks in one bar plot
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
