# ==========================================================================
# Script: 02_Cross_Species_Comparison_HFS_SHC.R
# Project: Wing Polyphenism in Hemipteran Insects
# Purpose: MetaNeighbor similarity, Harmony integration, and Top 50 Markers
# Species: Nilaparvata lugens (HFS) and Pyrrhocoris apterus (SHC)
# ==========================================================================

# 1. Load Required Libraries
library(SeuratDisk)
library(harmony)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(MetaNeighbor)
library(SingleCellExperiment)

# 2. Load Processed Data
# Ensure these files are in your 'data/' directory
load("data/hfs_obj_rename.Rdata")
hfs <- object
load("data/shc_obj_rename.Rdata")
shc <- object

# Set metadata for both species
hfs@meta.data$celltype <- Idents(hfs)
shc@meta.data$celltype <- Idents(shc)
hfs@meta.data$sp <- "NL" # Nilaparvata lugens
shc@meta.data$sp <- "PY" # Pyrrhocoris apterus

# 3. MetaNeighbor Analysis (Cell Type Similarity)
shc.sce <- as.SingleCellExperiment(shc)
hfs.sce <- as.SingleCellExperiment(hfs)

# Train model using HFS as reference
pretrained_model <- MetaNeighbor::trainModel(var_genes = rownames(hfs.sce),
                                             dat = hfs.sce,
                                             study_id = hfs.sce$sp,
                                             cell_type = hfs.sce$celltype)

# Run similarity test on SHC
aurocs <- MetaNeighborUS(trained_model = pretrained_model,
                         dat = shc.sce,
                         study_id = shc.sce$sp,
                         cell_type = shc.sce$celltype,
                         fast_version = TRUE)

# Save results to 'results/' directory
write.table(aurocs, quote = FALSE, "results/hfs_vs_shc.aurocs_celltype.xls", 
            sep = "\t", col.names = TRUE, row.names = TRUE)

pdf("results/hfs_vs_shc.similarity_heatmap.pdf", width = 7, height = 7)
plotHeatmapPretrained(aurocs, margins = c(10, 10))
dev.off()

# 4. Harmony Integration (Cross-species Batch Correction)
# Find common orthologous genes
common_genes <- intersect(rownames(shc), rownames(hfs))
shc_sub <- shc[common_genes, ]
hfs_sub <- hfs[common_genes, ]

# Merge datasets
pbmc <- merge(shc_sub, y = hfs_sub, add.cell.ids = c("PY", "NL"))
pbmc@meta.data$celltype1 <- paste(pbmc@meta.data$sp, pbmc@meta.data$celltype, sep = "_")

# Standard Seurat normalization and PCA
pbmc <- NormalizeData(pbmc) %>% 
        FindVariableFeatures(nfeatures = 2000) %>% 
        ScaleData() %>% 
        RunPCA(npcs = 30)

# Run Harmony to integrate species
pcSelect <- 20 # Adjust based on your ElbowPlot
pbmc <- RunHarmony(pbmc, group.by.vars = "sp", dims.use = 1:pcSelect)
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:pcSelect) %>%
        RunUMAP(reduction = "harmony", dims = 1:pcSelect)

# 5. Top 50 Marker Gene Visualization (Cross-species)
# Load the pre-calculated top 50 markers data
# Ensure 'sp.top50.Rdata' is in the 'data/' directory 
load("data/sp.top50.Rdata")

# Generate UMAP plots for Top 50 Markers analysis [cite: 10, 11]
# Note: 'cols' parameter is removed to use default R colors for compatibility
p_top1 <- DimPlot(pbmc, group.by = "celltype1", reduction = "umap", 
                  label = TRUE, repel = TRUE, label.box = TRUE) + 
          labs(title = "Top 50 Markers: Cell Type Clusters")
ggsave("results/TOP50.celltype.harmony.umap.pdf", p_top1, height = 15, width = 20, units = "cm")

p_top_split <- DimPlot(pbmc, group.by = "celltype1", reduction = "umap", 
                       label = TRUE, repel = TRUE, split.by = "sp") + 
               labs(title = "Top 50 Markers: Species Comparison")
ggsave("results/TOP50.celltype.harmony.split.umap.pdf", p_top_split, height = 15, width = 60, units = "cm")

p_top_sp <- DimPlot(pbmc, group.by = "sp", reduction = "umap", 
                    label = TRUE, repel = TRUE) + 
            labs(title = "Top 50 Markers: Species Distribution")
ggsave("results/TOP50.celltype.sp.harmony.umap.pdf", p_top_sp, height = 15, width = 20, units = "cm")
