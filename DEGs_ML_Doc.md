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
   - Code: **DEGxML_DESeq_analysis.R**
     - Input file: GSE226814_Mror_expr.txt (from: GEOarchive https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226814)
     - Output file: DEGs_biotrophic_vs_necrotrophic.csv

![PCA: Biotrophic vs. Necrotrophic](https://github.com/user-attachments/assets/2cd436b4-bfb7-4b91-8ea3-42b7c067abcd)
![VolcanoPlot: Biotrophic vs. Necrotrophic](https://github.com/user-attachments/assets/efe899ab-63e7-42da-9873-613cea1c79c0)
![Heatmap: Biotrophic vs. Necrotrophic](https://github.com/user-attachments/assets/4361b338-80e8-477a-9dc6-39ce5725d5d0)

3. Base dataset construction for Machine Learning
   - The database is constructed using the res_df obtained from the DESeq2 analysis. Rows = genes, Columns = features, DEG column = binary TAG for DEG
        - DEG binary tag criteria: 1 for padj < 0.05 & abs(log2FoldChange) >= 1, 0 for none.
        - Output file: DEGxML_dataset_base.csv

4. Extract additional features from gff3 and fasta files and merge the files in order to create a single file containing DEG stats and features.
   - Code: **DEGxML_gene_features_extraction.py** & **DEGxML_merge_datasets.py**
     - Input files: fasta -> MrorC26.groups.fasta, gff3 -> MrorC26.groups.gff3
     - Output file: DEGxML_features_from_gff_fasta.csv
   - Features extracted: gene_id, contig, start, end, strand, gene_len, n_exons, gc_content.
   - Merged files: DEGxML_features_from_gff_fasta.csv, DEGxML_dataset_base.csv
  
6. Data preprocessing for ML model and class balance verification
   - Code: **DEGxML_data_preprocessing.py** & **DEGxML_verify_class_balance.py**
     - Input file: DEGxML_merged_dataset.csv
     - Output files: DEGxML_features_scaled.csv, DEGxML_labels.csv
   - Exclude NaNs
   - Convert DEG to integers
   - Scaling values

Class distribution:
DEG
| 0 | 9369 |
| -- | -- |
| 1 | 1909 |

Class %:
DEG
| 0 | 83.1 % |
| -- | -- |
| 1 | 16.9 % |

7. Dataset splitting (80% train/20% test)
   - Code: **DEGxML_split.py**
     - Input files: DEGxML_features_scaled.csv, DEGxML_labels.csv
     - Output files: DEGxML_X_train.csv, DEGxML_X_test.csv, DEGxML_y_train.csv, DEGxML_y_test.csv 
