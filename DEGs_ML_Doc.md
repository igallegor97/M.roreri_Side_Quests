# DEGs prediction in *Moniliophthora roreri* isolates using Machine Learning

> For this side quest I'll be using publicly available RNA samples of *Moniliophthora roreri* that have been extracted from 2 different experimental conditions: 30dpi (biotrophic phase) & 60dpi (necrotrophic phase).
> 
> The data is stored at NCBI in the GEOarchive under the accession number GSE226814.

## Steps:

1. Data extraction
   - Download GSE226814_Mror_expr.txt.gz file and create the metadata table file **(DEGxML_metadata.tsv)**.
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

8. Model selection and training
   - Code: **DEGxML_model_training.py**
     - Input files: DEGxML_X_train.csv & DEGxML_y_train.csv
     - Output files: RandomForest_model.pkl, LogisticRegression_model.pkl, SVM_model.pkl, KNN_model.pkl
       
   - For this project 4 models will be tested:
     - Random Forest
     - Logistic Regression
     - SVM
     - KNN

       | Model | F1 score | std dev |
       | ---- | ---- | ---- |
       | Random Forest | 0.9993 | 0.0013 |
       | Logistic Regression | 0.9763 | 0.0062 |
       | SVM | 0.9609 | 0.0113 |
       | KNN | 0.8858 | 0.0066 |

10. Models testing and evaluation metrics
    - Code: **DEGxML_model_eval.py**
      - Input files: DEGxML_X_test.csv & DEGxML_y_test.csv
      - Output files: DEGxML_classification_reports.csv, roc_curve_KNN.png, roc_curve_LogisticRegression.png, roc_curve_RandomForest.png, roc_curve_SVM.png
        
    - The evaluation metrics are:
      - F1 score
      - Precision
      - Recall
      - Accuracy
      - Confusion matrix
      - AUC-ROC
      - 
| Model | Precision_No-DEG | Recall_No-DEG | F1_No-DEG | Precision_DEG | Recall_DEG | F1_DEG | Accuracy | AUC |
|--------------------|--------------------|--------------------|--------------------|--------------------|--------------------|--------------------|--------------------|--------------------|
| RandomForest       | 1.0                | 1.0                | 1.0                | 1.0                | 1.0                | 1.0                | 1.0                | 1.0                |
| LogisticRegression | 0.9983844911147012 | 0.9893276414087513 | 0.9938354328598231 | 0.949874686716792  | 0.9921465968586387 | 0.970550576184379  | 0.9898049645390071 | 0.9996186447780876 |
| SVM                | 0.9989183342347214 | 0.9855923159018143 | 0.9922105828632823 | 0.9336609336609336 | 0.9947643979057592 | 0.9632446134347274 | 0.987145390070922  | 0.9992344957450257 |
| KNN                | 0.9826933477555435 | 0.9695837780149413 | 0.9760945474080043 | 0.85995085995086   | 0.9162303664921466 | 0.8871989860583016 | 0.9605496453900709 | 0.9910144048902871 |

![confusion_matrix_RandomForest](https://github.com/user-attachments/assets/435db95e-5b02-4ec0-83df-805369b99809)
![confusion_matrix_LogisticRegression](https://github.com/user-attachments/assets/f48d683f-caba-4cf9-a7e2-39c3bbb34d59)
![confusion_matrix_SVM](https://github.com/user-attachments/assets/e13e1268-a706-4adf-b105-d62fedd6f881)
![confusion_matrix_KNN](https://github.com/user-attachments/assets/427dbf9b-00eb-4cba-b188-b792db9f4259)




