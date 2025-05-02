if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# Read file
counts <- read.delim("data/GSE226814_Mror_expr.txt", row.names = 1)

# Filter columns by condition selected 
selected_samples <- c("Mr30dpi_A", "Mr30dpi_B", "Mr30dpi_C",
                      "Mr60dpi_A", "Mr60dpi_B", "Mr60dpi_C")

counts_filtered <- counts[, selected_samples]

condition <- factor(c(rep("Biotrophic", 3), rep("Necrotrophic", 3)))

coldata <- data.frame(
  row.names = selected_samples,
  condition = condition
)

dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = coldata,
                              design = ~ condition)

# Filter low expression genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# DESeq2
dds <- DESeq(dds)

# Results
res <- results(dds, contrast = c("condition", "Necrotrophic", "Biotrophic"))

# Sort by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Save results
write.csv(as.data.frame(res_ordered), "results/DEGs_biotrophic_vs_necrotrophic.csv")

# PCA
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Calculate variance
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create plot with tags
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(title = "Biotrophic vs. Necrotrophic samples",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

# Create colum 'threshold' for gene classification
res_df$threshold <- "Not Significant"  
res_df$threshold[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$threshold[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

# Volcano plot

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.8, size = 2) +
  # Color: Up = red, Down = blue, No significance = gray
  scale_color_manual(values = c("Up" = "#E64A19", "Down" = "#1976D2", "Not Significant" = "grey70")) +
  # Cut lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  # Tags for most significant genes
  geom_text_repel(data = subset(res_df, padj < 1e-10 & abs(log2FoldChange) > 2),
                  aes(label = gene), size = 2.5, max.overlaps = 15) +
  # Title and tags
  labs(title = "Biotrophic vs Necrotrophic",
       x = "log2(Fold Change)",
       y = "-log10(padj)",
       color = "Expression") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Create conditions vector
conditions <- data.frame(
  Condition = factor(c("Biotrophic", "Biotrophic", "Biotrophic", 
                       "Necrotrophic", "Necrotrophic", "Necrotrophic"))
)
rownames(conditions) <- selected_samples

# Colors for each condition
ann_colors <- list(
  Condition = c(Biotrophic = "#1f78b4", Necrotrophic = "#e31a1c")
)

# Create heatmap with annotations and correlation numbers
pheatmap::pheatmap(cor_matrix,
                   annotation_col = conditions,
                   annotation_row = conditions,
                   annotation_colors = ann_colors,
                   display_numbers = TRUE,
                   number_format = "%.2f")

