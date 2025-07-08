# 🧬 DEG Prediction Using Genomic Features in *Moniliophthora roreri*

This repository contains a machine learning pipeline to predict whether a gene is differentially expressed (**DEG**) based solely on **genomic features**, avoiding expression-based bias or data leakage from RNA-Seq.

---

## 🎯 Project Goal

The goal is to predict DEG status using only gene-level genomic features, such as gene length, GC content, exon count, etc., derived from GFF3 annotations and genome sequences.

This pipeline avoids using expression values (e.g., TPM, counts) or DESeq2 results as input features, which would introduce **label leakage** into the model.

---

## 🛠️ Pipeline Overview

The workflow includes 3 main steps:

### 1. `feature_extraction_v2.py`

**Input:**
- A GFF3 file with gene annotations
- A genome FASTA file (optional, required for GC content)

**Output:**  
- `features/gene_features_v2.csv`

**Extracted features per gene:**
- `gene_length`
- `cds_length`
- `exon_count`
- `gc_content`
- `strand`
- `contig`
- `start`, `end`
- `start_position_norm` (normalized position within contig)
- `is_coding` (binary: has CDS or not)

---

### 2. `build_dataset_v2.py`

**Input:**
- `gene_features_v2.csv`
- DESeq2 output (`DEGxML_dataset_base.csv`)

**Processing:**
- Clean gene IDs (removing `.m01` suffix)
- Add DEG label:
  - `padj < 0.05`
  - `|log2FoldChange| >= 1`
- Merge features with labels by `gene_id`

**Output:**
- `features/deg_dataset_v2.csv`

This file contains genomic features and a binary `DEG` column (`1 = DEG`, `0 = not DEG`).

---

### 3. `train_RF_v2.py`

**Purpose:** Train and evaluate a Random Forest classifier using only genomic features.

**Steps:**
- Load dataset
- Remove non-numeric columns (e.g., `contig`, `strand`)
- Fill missing values (e.g., `NaN` in `gc_content`) with column means
- Train/test split (`80/20`, stratified)
- Fit a `RandomForestClassifier` with class weighting
- Evaluate performance:
  - Classification report (precision, recall, F1)
  - Confusion matrix
  - ROC AUC score
- Interpretability:
  - Feature importances
  - SHAP value analysis

**Outputs:**
- `features/feature_importance_RF.png`
- `features/shap_summary_RF.png`
- Console logs with metrics and evaluation results

---

## 🧪 Biological Rationale

This approach allows us to explore whether **genomic characteristics** are associated with differential expression, without using RNA-Seq values as input.

It is particularly useful for:
- Transfer learning between strains or conditions
- Identifying genomic predictors of regulation
- Understanding non-expression-based regulation patterns

---

## 📁 Folder Structure

M.roreri_Side_Quests/
│
├── data/
│ └── DEGxML_dataset_base.csv # DESeq2 output
│
├── features/
│ ├── gene_features_v2.csv # Extracted genomic features
│ ├── deg_dataset_v2.csv # Features + DEG label
│ ├── feature_importance_RF.png # Feature importance plot
│ └── shap_summary_RF.png # SHAP summary plot
│
├── feature_extraction_v2.py
├── build_dataset_v2.py
├── train_RF_v2.py
└── README.md


---

## 🔮 Future Directions

- Use other ML models (XGBoost, Logistic Regression, SVM)
- Perform hyperparameter tuning and cross-validation
- Extend to multi-condition DEG classification
- Build models to predict **pangenome category** (core, accessory, unique)
- Integrate functional annotation or synteny data

---

## 📌 Dependencies

Make sure to install:

```bash
pip install pandas scikit-learn matplotlib seaborn shap

👩‍🔬 Contact
Project developed by Isabella Gallego Rendón
Master’s Student in Bioinformatics – University of São Paulo 🇧🇷
RSG Colombia Chair – ISCB Student Council
GitHub: @igallegor97
