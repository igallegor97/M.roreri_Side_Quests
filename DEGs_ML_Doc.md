# DEGs prediction in *Moniliophthora roreri* isolates using Machine Learning

> For this side quest I'll be using publicly available RNA samples of *Moniliophthora roreri* that have been extracted from 2 different experimental conditions: 30dpi (biotrophic phase) & 60dpi (necrotrophic phase).
> 
> The data is stored at NCBI in the GEOarchive under the accession number GSE226814.

## Steps:

1. Download GSE226814_Mror_expr.txt.gz file and create the metadata table file **(DEGxML_metadata.tsv)**.
   - Table structure preview:

| sample_ID | condition |
| --------- | --------- |
| Mr30dpi_A | biotrophic |
| Mr30dpi_B | biotrophic |
| Mr30dpi_C | biotrophic |
| Mr60dpi_A | necrotrophic |
| Mr60dpi_B | necrotrophic |
| Mr60dpi_C | nectrotrophic |

2. DESeq2 analysis
   - Code: DEGxML_DESeq_analysis.R

![PCA: Biotrophic vs. Necrotrophic](https://github.com/user-attachments/assets/2cd436b4-bfb7-4b91-8ea3-42b7c067abcd)
![VolcanoPlot: Biotrophic vs. Necrotrophic](https://github.com/user-attachments/assets/efe899ab-63e7-42da-9873-613cea1c79c0)
![Heatmap: Biotrophic vs. Necrotrophic](https://github.com/user-attachments/assets/4361b338-80e8-477a-9dc6-39ce5725d5d0)

3. Base dataset construction for Machine Learning
   - The database is constructed using the res_df obtained from the DESeq2 analysis. Rows = genes, Columns = features, DEG column = binary TAG for DEG
        - DEG binary tag criteria: 1 for padj < 0.05 & abs(log2FoldChange) >= 1, 0 for none.
