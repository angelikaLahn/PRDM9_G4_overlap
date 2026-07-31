# ==============================================================
# Canonical G4 motif detection 
# ==============================================================
# Author: Dr. Angelika Lahnsteiner
# Date: 2026-07-30
# ==============================================================

# ------------------------------
# Load libraries
# ------------------------------
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Biostrings)
library(GenomicRanges)
library(data.table)


# Define file path
pqsfinder_hg38 <- read.delim("~/path_to_predicted_g4s", header=FALSE)

# mutate all to upper case
pqsfinder_hg38 <- pqsfinder_hg38 %>%
  mutate(V7 = toupper(as.character(V7)))

# ==============================================================
# Reverse-complement sequences if C-rich (C > G)
# ==============================================================
# Prepare just the G-motifs, reverse complement in the case of C-rich regions
# Count C and G in all sequences at once
rc_if_C_rich <- function(seq_vec) {
  c_count <- str_count(seq_vec, "C")
  g_count <- str_count(seq_vec, "G")
  to_rc <- c_count > g_count
  seq_vec[to_rc] <- as.character(reverseComplement(DNAStringSet(seq_vec[to_rc])))
  return(seq_vec)
}

bed_data <- pqsfinder_hg38 %>%
  mutate(
    sequence_comp  = rc_if_C_rich(V7),
  )


bed_data <- bed_data[c(1,2,3,5,8)]
colnames(bed_data)[4] <- "sequence"
colnames(bed_data)

# View first few rows
head(bed_data)

# ==============================================================
# 1. Define canonical G4 pattern (strict canonical)
# ==============================================================
# Four G-runs (≥3 Gs each), three loops (1–12 nt), 
# loops may include Gs or GG but not GGG tracts.
canonical_pattern_short <- paste0(
  "(G{3,4})",                # stem1
  "([ATCG]{1,3})",           # loop1
  "(G{3,4})",                # stem2
  "([ATCG]{1,3})",           # loop2
  "(G{3,4})",                # stem3
  "([ATCG]{1,3})",           # loop3
  "(G{3,4})"                # stem4
)

#filter
G4s_short_loops <- bed_data[str_detect(bed_data$sequence_comp, canonical_pattern_short), ]

canonical_pattern_medium <- paste0(
  "(G{3,4})",                # stem1
  "([ATCG]{4,6})",           # loop1
  "(G{3,4})",                # stem2
  "([ATCG]{4,6})",           # loop2
  "(G{3,4})",                # stem3
  "([ATCG]{4,6})",           # loop3
  "(G{3,4})"                # stem4
)

#filter
G4s_medium_loops <- bed_data[str_detect(bed_data$sequence_comp, canonical_pattern_medium), ]


canonical_pattern_long <- paste0(
  "(G{3,4})",                # stem1
  "([ATCG]{7,12})",           # loop1
  "(G{3,4})",                # stem2
  "([ATCG]{7,12})",           # loop2
  "(G{3,4})",                # stem3
  "([ATCG]{7,12})",           # loop3
  "(G{3,4})"                # stem4
)

#filter
G4s_long_loops <- bed_data[str_detect(bed_data$sequence_comp, canonical_pattern_long), ]



# check if there are any motifs shared
all_motifs <- c(
  G4s_short_loops$sequence_comp,
  G4s_medium_loops$sequence_comp,
  G4s_long_loops$sequence_comp
)

motif_counts <- table(all_motifs)

#extract the motifs which are shared in at least two datasets
shared_motifs <- names(motif_counts[motif_counts >= 2])
#prepare a dataframe
shared_motifs_df <- data.table(sequence_comp = common_motifs)

# remove the shared once from the original tables
G4s_short_loops <- G4s_short_loops %>% filter(!sequence_comp %in% shared_motifs,)
G4s_medium_loops <- G4s_medium_loops %>% filter(!sequence_comp %in% shared_motifs,)
G4s_long_loops <- G4s_long_loops %>% filter(!sequence_comp %in% shared_motifs,)

colnames(G4s_short_loops) <- c("chr", "start", "end", "motif_length", "G4_motif_sequence")
colnames(G4s_medium_loops) <- c("chr", "start", "end", "motif_length", "G4_motif_sequence")
colnames(G4s_long_loops) <- c("chr", "start", "end", "motif_length", "G4_motif_sequence")


write.table(G4s_short_loops, "~/Desktop/G4_PRDM9_Analysis/G4s/Canonical_G4s/G4s_short_loops.tsv", quote=FALSE, row.names = FALSE, sep=",")
write.table(G4s_medium_loops, "~/Desktop/G4_PRDM9_Analysis/G4s/Canonical_G4s/G4s_medium_loops.tsv", quote=FALSE, row.names = FALSE, sep=",")
write.table(G4s_long_loops, "~/Desktop/G4_PRDM9_Analysis/G4s/Canonical_G4s/G4s_long_loops.tsv", quote=FALSE, row.names = FALSE, sep=",")
write.table(shared_motifs_df, "~/Desktop/G4_PRDM9_Analysis/G4s/Canonical_G4s/G4s_nested_canonical_motifs.tsv", quote=FALSE, row.names = FALSE, sep=",")
