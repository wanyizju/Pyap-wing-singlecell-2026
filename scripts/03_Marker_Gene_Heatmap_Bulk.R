# ==========================================================================
# Script: 03_Marker_Gene_Heatmap_Bulk.R
# Purpose: Visualization of Marker Genes in Bulk RNA-seq
# ==========================================================================

library(pheatmap)
library(dplyr)
library(readr)

# 1. Load Data
bulk_expr <- read_csv("data/bulk_expression_complete.csv")
marker_df <- read_delim("data/marker_gene_list.txt", delim = "\t")

# 2. Filter and Transform
rownames(bulk_expr) <- bulk_expr$GeneID
log_expr <- bulk_expr[bulk_expr$GeneID %in% marker_df$GeneID, ] %>%
            mutate(InR2_mean = rowMeans(select(., py_InR2_1, py_InR2_2, py_InR2_3)),
                   WT_mean = rowMeans(select(., py_wt_1, py_wt_2, py_wt_3))) %>%
            select(InR2_mean, WT_mean) %>%
            log2() %>%
            as.matrix()

# 3. Annotation Setup
marker_anno <- marker_df[match(rownames(log_expr), marker_df$GeneID), ]
rownames(log_expr) <- marker_anno$GeneSymbol
annotation_row <- data.frame(CellType = marker_anno$CellType)
rownames(annotation_row) <- marker_anno$GeneSymbol

# 4. Plot Heatmap
pheatmap(log_expr, 
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = annotation_row,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "results/Bulk_Marker_Heatmap.pdf", width = 6, height = 10)
