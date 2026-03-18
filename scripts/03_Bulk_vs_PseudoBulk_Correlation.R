# ==========================================================================
# Script: 03_Bulk_vs_PseudoBulk_Correlation.R
# Purpose: Correlation between Bulk RNA-seq and Pseudo-bulk from scRNA-seq
# ==========================================================================

library(ggplot2)
library(dplyr)
library(readr)

# 1. Load Data (Using relative paths)
bulk_expr <- read_csv("data/bulk_expression_complete.csv")
pseudo_expr <- read_csv("data/pseudo_bulk_expression.csv")

# 2. Processing Bulk Data
# Calculate mean for the experimental group (py_InR2)
bulk_expr <- bulk_expr %>%
  mutate(Exp_mean = rowMeans(select(., py_InR2_1, py_InR2_2, py_InR2_3)))

# 3. Processing Pseudo-bulk Data
common_genes <- intersect(bulk_expr$GeneID, pseudo_expr$GeneID)
bulk_common <- bulk_expr %>% filter(GeneID %in% common_genes)
pseudo_common <- pseudo_expr %>% filter(GeneID %in% common_genes)

pseudo_common_numeric <- pseudo_common %>% select(where(is.numeric))
pseudo_common$Pseudo_mean <- rowMeans(pseudo_common_numeric)

# 4. Merge and Log2 Transformation
plot_data <- merge(bulk_common[, c("GeneID", "Exp_mean")],
                   pseudo_common[, c("GeneID", "Pseudo_mean")], by = "GeneID") %>%
             mutate(log_Bulk = log2(Exp_mean + 1),
                    log_Pseudo = log2(Pseudo_mean + 1))

# 5. Pearson Correlation and Plotting
r_val <- round(cor(plot_data$log_Pseudo, plot_data$log_Bulk, method = "pearson"), 2)

p <- ggplot(plot_data, aes(x = log_Pseudo, y = log_Bulk)) +
  geom_point(alpha = 0.4, size = 1.5, color = "purple") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(x = "log2(Pseudo-bulk Expression + 1)",
       y = "log2(Bulk Expression + 1)",
       title = paste0("Bulk vs Pseudo-bulk Correlation (r = ", r_val, ")")) +
  theme_minimal()

ggsave("results/Bulk_PseudoBulk_Correlation.pdf", p, width = 6, height = 5)