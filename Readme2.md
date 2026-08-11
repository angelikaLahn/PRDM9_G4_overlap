# Analysis of G-quadruplex DNA structures at PRDM9 binding sites

## Overview

This repository contains the custom scripts used for the computational analyses described in:

**Hartge MK, Hernandez NP, Pavlik T, Polakova N, Schoenauer E, Huber J, Striedner Y, Mair T, Tiemann-Boege I, Weissensteiner MH, Brandstetter J, Risch A, Cechova M, Lahnsteiner A.**  
**G-quadruplex structures act as a novel recognition motif for the meiosis-specific histone methyltransferase PRDM9.**  
bioRxiv (2026). DOI: `10.64898/2026.07.31.742039`

The analyses investigate the occurrence of experimentally detected and predicted **G-quadruplex (G4) DNA structures at PRDM9 binding sites** in human HEK293 cells.

In particular, the scripts were used to:

1. identify reproducible G4 CUT\&Tag peaks across datasets;
2. annotate G4 and PRDM9 peaks to genomic features;
3. investigate the strand orientation of G4 motifs relative to PRDM9 binding sites;
4. identify genomic regions shared between G4 CUT\&Tag and PRDM9 ChIP-seq peaks; and
5. characterize canonical G4 motifs according to loop length.

\---

## Repository structure

```text
PRDM9_G4_overlap/
│
├── datasets/
│   ├── SEACR_output_G4s/
│   │   └── G4 CUT&Tag peak files from individual datasets
│   │
│   └── ChIPseeker/
│       ├── HEK293_G4s_hg38.bed
│       ├── HEK293_PRDM9-B_HA_MACS2_hg38.bed
│       └── G4_PRDM9_25bp_window.bed
│
├── scripts/
│   ├── 01_identify_shared_G4_peaks.sh
│   ├── 02_ChIPseeker_annotation.R
│   ├── 03_Analyse_g4_strand_orientation.R
│   ├── 04_G4_PRDM9_window_intersect.sh
│   └── 05_canonical_G4_loop_lengths.R
│
├── LICENSE
└── README.md
```

\---

# Analysis workflow

## 1\. Peak calling

Peak calling of the G4 CUT\&Tag and PRDM9 ChIP-seq datasets was performed using the **Galaxy platform**.

### G4 CUT\&Tag

G4 CUT\&Tag peaks were called using **SEACR**.

Peak files from the individual datasets are provided in:

```text
datasets/SEACR_output_G4s/
```

The original SRA accession numbers are retained in the filenames.

### PRDM9 ChIP-seq

PRDM9 ChIP-seq peaks were called using **MACS2**, using the corresponding input dataset as background.

The resulting PRDM9 peak file used for the downstream analyses is:

```text
datasets/ChIPseeker/HEK293_PRDM9-B_HA_MACS2_hg38.bed
```

All genomic coordinates used in this repository correspond to the **GRCh38/hg38 human reference genome**.

\---

## 2\. Identification of reproducible G4 peaks

To generate a consensus set of experimentally supported G4 regions, peaks detected in the individual G4 CUT\&Tag datasets were compared using:

```bash
bedtools multiinter
```

The analysis is described in:

```text
scripts/01_identify_shared_G4_peaks.sh
```

Regions detected in **at least two G4 CUT\&Tag datasets** were retained.

The resulting consensus G4 peak file is:

```text
datasets/ChIPseeker/HEK293_G4s_hg38.bed
```

\---

## 3\. Annotation of G4 and PRDM9 peaks to genomic features

G4 and PRDM9 peaks were annotated to genomic features using the **ChIPseeker** R/Bioconductor package.

The analysis is implemented in:

```text
scripts/02_ChIPseeker_annotation.R
```

The following peak files are used as input:

```text
HEK293_G4s_hg38.bed
HEK293_PRDM9-B_HA_MACS2_hg38.bed
```

Peak distributions are annotated relative to genomic features including:

* promoters / transcription start sites;
* exons;
* introns;
* downstream regions; and
* intergenic regions.

Transcription start site regions were defined as ±3 kb around the TSS.

\---

## 4\. G4 strand-orientation analysis

To investigate whether PRDM9 binding shows a preference for G4 motifs located on either DNA strand, experimentally detected G4 regions were assigned to nearby predicted G4 motifs.

The analysis is provided in:

```text
scripts/03_Analyse_g4_strand_orientation.R
```

For each G4 CUT\&Tag peak, the nearest predicted G4 motif was identified using `GenomicRanges`.

Only predicted G4 motifs located within **±500 bp of the G4 peak center** were retained.

The strand information of the corresponding predicted G4 motif was subsequently used to investigate PRDM9 binding relative to G4 orientation.

### Additional input

This analysis requires a BED file containing the genomic positions, strand information, and sequence of predicted G4 motifs generated using **pqsfinder**.

The file can be directly downloaded from the Website (https://pqsfinder.fi.muni.cz/genomes) of pqsfinder with hg38 annotation, or generated in R using the pqsfinder package. The path to this file needs to be specified in the script before running the analysis.

\---

## 5\. Identification of shared G4 and PRDM9 regions

To determine whether experimentally detected G4 structures occur in close proximity to PRDM9 binding sites, G4 CUT\&Tag peaks were compared with PRDM9 ChIP-seq peaks using:

```bash
bedtools window
```

The commands are provided in:

```text
scripts/04_G4_PRDM9_window_intersect.sh
```

A PRDM9 peak was considered associated with a G4 peak when it occurred within a **25-bp window upstream or downstream of the G4 peak**:

```bash
bedtools window \\
    -a HEK293_G4s_hg38.bed \\
    -b HEK293_PRDM9-B_HA_MACS2_hg38.bed \\
    -w 25 \\
    > G4_PRDM9_25bp_window.bed
```

The resulting file is provided as:

```text
datasets/ChIPseeker/G4_PRDM9_25bp_window.bed
```

\---

## 6\. Classification of canonical G4 motifs by loop length

Canonical G4 motifs were further classified according to their loop lengths using:

```text
scripts/05_canonical_G4_loop_lengths.R
```

Predicted G4 sequences were first oriented according to their G-rich strand. C-rich sequences were reverse-complemented before motif classification.

Canonical G4s were defined as sequences containing four G-runs with 3–4 guanines per run.

G4s were subsequently divided into three loop-length classes:

|Class|Loop length|
|-|-:|
|Short loops|1–3 nt|
|Medium loops|4–6 nt|
|Long loops|7–12 nt|

This analysis was used to investigate whether the association between predicted G4 motifs and PRDM9 binding sites depends on G4 loop length.

\---

# Dependencies

The analyses were performed using **R ≥ 4.5.0** and **bedtools ≥ 2.30.0**.

### R / Bioconductor packages

The scripts use the following packages:

```text
ChIPseeker
TxDb.Hsapiens.UCSC.hg38.knownGene
org.Hs.eg.db
clusterProfiler
GenomicRanges
rtracklayer
Biostrings
stringr
dplyr
tidyr
data.table
ggplot2
ggpubr
```

Predicted G4 motifs were generated using **pqsfinder**.

### Reference genome

```text
GRCh38 / hg38
```

\---

# Running the analyses

Clone the repository:

```bash
git clone https://github.com/angelikaLahn/PRDM9_G4_overlap.git
cd PRDM9_G4_overlap
```

The scripts are numbered according to the approximate order of the analyses:

```text
01_identify_shared_G4_peaks.sh
        ↓
02_ChIPseeker_annotation.R
        ↓
03_Analyse_g4_strand_orientation.R
        ↓
04_G4_PRDM9_window_intersect.sh
        ↓
05_canonical_G4_loop_lengths.R
```

R scripts can be executed, after adjusting the required input paths, using for example:

```bash
Rscript scripts/02_ChIPseeker_annotation.R
```

and:

```bash
Rscript scripts/03_Analyse_g4_strand_orientation.R
```

The `.sh` files contain the corresponding `bedtools` command-line workflows.

**Note:** Some scripts currently contain local example paths. These need to be replaced by paths corresponding to the location of the cloned repository before execution.

\---

# Data availability

Processed BED files required for the overlap analyses are included in the `datasets/` directory.

The individual G4 CUT\&Tag datasets are identified by their original **SRA accession numbers** in the filenames under:

```text
datasets/SEACR_output_G4s/
```

Raw sequencing data can be obtained from the respective public sequence repositories referenced in the manuscript.

\---

# Citation

If you use this repository or the associated analysis workflow, please cite:

> Hartge MK, Hernandez NP, Pavlik T, Polakova N, Schoenauer E, Huber J, Striedner Y, Mair T, Tiemann-Boege I, Weissensteiner MH, Brandstetter J, Risch A, Cechova M, Lahnsteiner A.  
> \*\*G-quadruplex structures act as a novel recognition motif for the meiosis-specific histone methyltransferase PRDM9.\*\*  
> bioRxiv (2026). DOI: `10.64898/2026.07.31.742039`

\---

# License

This project is licensed under the **MIT License**. See the `LICENSE` file for details.

\---

# Contact

**Dr. Angelika Lahnsteiner**  
Department of Biosciences and Medical Biology  
Paris Lodron University of Salzburg  
Salzburg, Austria

