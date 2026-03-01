# Pyap-wing-singlecell-2026

This repository contains the R scripts and processed data for Transcription Factor (TF) analysis in *Pyrrhocoris apterus* (SHC) and *Nilaparvata lugens* (HFS).

## 1. System Requirements
- **R Version:** >= 4.2.0
- **Operating System:** Windows, macOS, or Linux.
- **Required Packages:** `SCENIC`, `ComplexHeatmap`, `SCopeLoomR`, `circlize`, `data.table`, `dplyr`.

## 2. Repository Structure
- `scripts/`: 
  - `05_TF_Regulon_Activity_Visualization.R`: Main script for visualizing TF activity.
- `data/`: 
  - `gy_vs_shc.homogene.xls`: Gene homology mapping table for *P. apterus*.
  - `gy_vs_hfs.homogene.xls`: Gene homology mapping table for *N. lugens*.
  - `SHC_SCENIC.loom`: Regulon AUC data (processed by pySCENIC).

## 3. Installation Guide
To install the necessary R packages, run:
```R
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SCENIC", "ComplexHeatmap", "SCopeLoomR"))
