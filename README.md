# Analysis of G-quadruplex DNA structure occurrence at PRDM9 binding sites in HEK293 cells


Overview

This repository contains scripts used to identify shared regions of G-quadruplex (G4) CUT&Tag peaks with PRDM9 ChIP-seq peaks to assess their concurrence. The following steps have been performed: 

1. Peak calling
Peak calling of BG4 CUT&Tag and PRDM9 ChIP-seq datasets was performed on the Galaxy platform. In short, xx
The following datasets have been used:

* G4s:
* PRDM9: 

The individual peak datasets from SEACR for BG4 CUT&Taq are saved in the folder sample_files/SEACR_output/.
The PRDM9 peak file from two replicates compared to a input background file in MACS2 was saved in sample_files/HEK293_PRDM9-B_HA_MACS2_hg38.bed.


2. Identify BG4 CUT&Tag peaks shared of at least two replicates

Using bedtools multiinter, we identified peaks which are shared by at least two replicates saved as bedfiles in sample_files/SEACR_output/ using the script 01_identify_shared_G4_peaks.md. The resulting file is stored as sample_files/HEK293_G4s_hg38.bed. 

3. Annotation to genomic features using ChIP-seeker R-package

To annotate the resulting peak files HEK293_G4s_hg38.bed for G4s and HEK293_PRDM9-B_HA_MACS2_hg38.bed for PRDM9 to genomic features such as introns, exons, intergenic regions or the transcription start site, we have used the ChIPseeker package in R using the custom script 02_ChIPseeker_annotation.R 
To analyse a potential strand bias, we plotted the occurence of PRDM9 peaks on the positive or negative strand using the G4 motifs as centers (03_Analyse_g4_strand_orientation.R). 

4. Identify shared peaks from BG4 CUT&Tag and PRDM9 ChIP-seq datasets

The overlap of PRDM9 binding sites with G4 forming sites was assessed using bedtools window with a 25-bp window up- and downstream of the G4 peak.  

5. Predict canonical G4s with different loop lengths
To predict canonical G4s meeting the expression 

7. Dependencies

R ≥ 4.5.0
pqsfinder
ChIPseeker
GenomicRanges
annotatr
clusterProfiler
org.Hs.eg.db
bedtools ≥ 2.30.0
Reference genome: GRCh38

Sources:

