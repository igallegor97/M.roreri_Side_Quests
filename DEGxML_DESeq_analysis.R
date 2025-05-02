if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")

library(DESeq2)
library(tidyverse)

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
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

# Volcano plot
res_df <- as.data.frame(res_ordered)
res_df$gene <- rownames(res_df)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed") +
  labs(title = "Volcano plot", x = "log2 Fold Change", y = "-log10 Adjusted P-Value")
