# Pyap-wing-singlecell-2026

This repository contains the R scripts and processed data for the study of insect wing development in *Pyrrhocoris apterus* (SHC) and *Nilaparvata lugens* (HFS).

## 1. System Requirements
- **R Version:** >= 4.2.0
- **Operating System:** Windows, macOS, or Linux.
- **Key R Packages:** `Seurat`, `harmony`, `MetaNeighbor`, `SCENIC`, `ComplexHeatmap`, `pheatmap`, `tidyverse`.

## 2. Repository Structure
### `scripts/`
- **`01_Cross_Species_Comparison_HFS_SHC.R`**: Performs cell type similarity assessment (MetaNeighbor) and cross-species integration (Harmony), including Top 50 marker gene analysis.
- **`02_Bulk_vs_PseudoBulk_Correlation.R`**: Calculates the Pearson correlation between scRNA-seq pseudo-bulk data and traditional Bulk RNA-seq data.
- **`03_Marker_Gene_Heatmap_Bulk.R`**: Visualizes the expression of scRNA-seq defined marker genes within Bulk RNA-seq datasets.
- **`04_TF_Regulon_Activity_Visualization.R`**: Visualizes Transcription Factor (TF) regulon activity scores derived from pySCENIC.

### `data/`
- **Expression Matrices**: `pseudo_bulk_expression.csv`.
- **Mapping & Lists**: `gy_vs_shc.homogene.xls`, `marker_gene_list.txt`.

### `results/`
- All scripts are configured to automatically save output figures (PDF/PNG) and tables (XLS) into this directory.

## 3. Data Availability (Large Files)
Due to GitHub's file size limitations, the following core Seurat objects are hosted on **Zenodo**:
- **Files**: `hfs_obj_rename.Rdata`, `shc_obj_rename.Rdata`.
- **DOI/Link**: [Insert your Zenodo DOI Link here]
- **Instruction**: Please download these files and place them into the `data/` folder before running script `01`.

## 4. Usage Instructions
1. **Setup**: Clone the repository and ensure the `results/` folder exists.
2. **Data**: Place all required data files (including those from Zenodo) into the `data/` directory.
3. **Execution**: Run the scripts in numerical order (`01` to `04`) to reproduce the findings presented in the manuscript.

**Estimated Run Time**: Each script typically completes in < 10 minutes on a standard computational workstation.
