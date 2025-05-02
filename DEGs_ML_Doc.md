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
