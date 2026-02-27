# ==========================================================================
# Script: 05_TF_Regulon_Activity_Visualization.R
# Project: Wing Polyphenism in Pyrrhocoris apterus
# Purpose: Visualization of Transcription Factor (TF) Regulon Activity
# ==========================================================================

# 1. Load Required Libraries
library(SCopeLoomR)
library(SCENIC)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(dplyr)

# 2. Environment and Data Setup
# Use relative paths for GitHub compatibility. 
# Ensure these files are placed in a 'data/' directory within your repository.
scenicLoomPath <- "data/SHC_SCENIC.loom" 
[cite_start]homologyTablePath <- "data/gy_vs_shc.homogene.xls" [cite: 1]

# Load the homology mapping table
# This file must contain 'gene_id' and 'ID2' columns
[cite_start]rt1 <- read.table(homologyTablePath, sep = "\t", header = TRUE) [cite: 3]

# Open loom file and extract Regulon AUC data
[cite_start]loom <- open_loom(scenicLoomPath) [cite: 1]
[cite_start]regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC') [cite: 1]
[cite_start]close_loom(loom) [cite: 1]

# 3. Data Processing
# Keep only non-duplicated extended regulons
[cite_start]sub_regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),] [cite: 2]

# Calculate average activity per cell type
# [cite_start]'cellsPerGroup' should be defined based on your Seurat/SingleCellExperiment clusters [cite: 1, 2]
regulonActivity_byGroup <- sapply(cellsPerGroup, function(cells) 
  [cite_start]rowMeans(getAUC(sub_regulonAUC)[,cells])) [cite: 2]

# Scale and clean data
[cite_start]df_scaled <- data.frame(t(scale(t(regulonActivity_byGroup), center = TRUE, scale = TRUE))) [cite: 3]
[cite_start]df_scaled <- na.omit(df_scaled) [cite: 3]

# [cite_start]Rename and reorder columns for publication consistency [cite: 3]
colnames(df_scaled) <- c("Epithelial cell", "Tracheal cell", "Neuron cell", 
                         "Hemocyte", "Glial cell", "Muscle cell")
[cite_start]df_scaled <- df_scaled[, c(1, 2, 3, 5, 4, 6)] # Reorder based on your study design [cite: 3]

# 4. Map Homology Information
# [cite_start]Clean TF names and merge with homology IDs [cite: 3]
[cite_start]df_scaled$ID2 <- gsub("\\(\\+\\)", "", rownames(df_scaled)) [cite: 3]
[cite_start]rt4 <- merge(df_scaled, rt1, by = "ID2", sort = FALSE) [cite: 3]
[cite_start]rt4$ID4 <- paste(rt4$ID2, " (", rt4$gene_id, ")", sep = "") [cite: 3]
rownames(rt4) <- rt4$ID4
[cite_start]rt5 <- rt4[, 1:6] # Retain only numeric columns for the heatmap [cite: 3]

# 5. Generate Heatmap Output
[cite_start]color_palette <- c("#9777cc", "white", "#d15047") [cite: 3]

# Save as PDF for high-quality publication figures
[cite_start]pdf("results/TF_Activity_Heatmap.pdf", height = 10, width = 7) [cite: 3, 8]
Heatmap(as.matrix(rt5), 
        [cite_start]name = "TF regulon activity", [cite: 4]
        [cite_start]col = color_palette, [cite: 4]
        show_row_names = TRUE, 
        [cite_start]show_column_names = TRUE, [cite: 4, 5]
        [cite_start]row_names_gp = gpar(fontsize = 8), [cite: 9]
        [cite_start]column_names_gp = gpar(fontsize = 8), [cite: 9]
        [cite_start]clustering_method_rows = "ward.D2", [cite: 4]
        [cite_start]cluster_columns = FALSE, [cite: 4, 10]
        [cite_start]column_names_rot = 45) [cite: 4, 12]
dev.off()