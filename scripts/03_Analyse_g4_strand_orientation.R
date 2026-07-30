# ==============================================================================
#      Analysis of strand specific binding of PRDM9 to G4 structures
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
library(rtracklayer)
library(GenomicRanges)

# ------------------------------
# Set directories and load files
# ------------------------------
setwd("~/Desktop/G4_PRDM9_Analysis/Analysis/ChIPseeker")
samplefiles <- list.files("~/Desktop/G4_PRDM9_Analysis/Analysis/ChIPseeker/files", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
print(samplefiles)
names(samplefiles) <- c("G4s","PRDM9","pqs")
print(samplefiles)

#import as peak file
pqs <- rtracklayer::import(samplefiles[["pqs"]])
g4 <- rtracklayer::import(samplefiles[["G4s"]])
prdm9 <- rtracklayer::import(samplefiles[["PRDM9"]])

#import as bed files
#g4 <- read.delim("~/Desktop/G4_PRDM9_Analysis/Analysis/ChIPseeker/files/HEK293_G4s_hg38.bed", header=FALSE)
#prdm9 <- read.delim("~/Desktop/G4_PRDM9_Analysis/Analysis/ChIPseeker/files/HEK293_PRDM9-B_HA_MACS2_hg38.bed", header=FALSE)
#pqs <- read.delim("~/Desktop/G4_PRDM9_Analysis/Analysis/ChIPseeker/files/pqsfinder_hg38.bed", header=FALSE)

#check length
length(g4)
length(pqs)
length(prdm9)

table(strand(pqs))

# ------------------------------
# Analysis
# ------------------------------
# find the center of the G4 peaks
g4_center <- resize(g4, width = 1, fix = "center")

hits <- nearest(g4_center, pqs, ignore.strand = TRUE)
dist <- distance(g4_center, pqs[hits])
# NAs detected --> peaks contain no G4 motif --> only 267 peaks in total
sum(is.na(hits))
[1] 267
> length(hits)
[1] 43557

valid <- !is.na(hits)

#create a mask with all valid entries
dist <- distance(g4_center, pqs[hits])
g4_valid <- g4_center[valid]
hits_valid <- hits[valid]
pqs_valid <- pqs[hits_valid]

dist <- distance(g4_valid, pqs_valid)

# keep all G4s with a distance of +/-500bp to the peak center
keep <- which(dist <= 500)

g4_filtered <- g4_valid[keep]
pqs_matched <- pqs_valid[keep]

strand(g4_filtered) <- strand(pqs_matched)

g4_plus  <- g4_filtered[strand(g4_filtered) == "+"]
g4_minus <- g4_filtered[strand(g4_filtered) == "-"]

# Now assign the prdm9 peaks to the G4_plus and G4_minus datasets
g4_plus_center  <- resize(g4_plus,  width = 1, fix = "center")
g4_minus_center <- resize(g4_minus, width = 1, fix = "center")

#build symmetric windows around the centers
# 3. Extend the windows upstream and downstream (e.g., +/- 3000 bp)
custom_windows_plus <- makeBioRegionFromGranges(
  gr = g4_plus, 
  by = "g4_plus",            
  type = "start_site", 
  upstream = 5000, 
  downstream = 5000
)

custom_windows_minus <- makeBioRegionFromGranges(
  gr = g4_minus, 
  by = "g4_minus",            
  type = "start_site", 
  upstream = 5000, 
  downstream = 5000
)

# 4. Generate the matrix using your custom windows
tagMatrix_plus <- getTagMatrix(prdm9, windows = custom_windows_plus)
tagMatrix_minus <- getTagMatrix(prdm9, windows = custom_windows_minus)
# 5. Plot the profile centered on G4 locations
plotAvgProf(tagMatrix_plus, xlim = c(-5000, 5000), xlab = "G4 Center", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix_minus, xlim = c(-5000, 5000), xlab = "G4 Center", ylab = "Read Count Frequency")


# ------------------------------
# Plot the heatmap
# ------------------------------

# prepare a custom heatmap
p <- tagHeatmap(tagMatrix_plus, xlab = "G4 Center")
# add red and white for plus strand
p + scale_fill_gradient(low = "white", high = "pink") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )


# prepare a custom heatmap
p1 <- tagHeatmap(tagMatrix_minus, xlab = "G4 Center")
# add blue and white for minus strand
p1 + scale_fill_gradient(low = "white", high = "lightblue") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

