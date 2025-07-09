# Promoter-Based DEG Classification Using CNN

## Overview

This repository contains the workflow and scripts for classifying Differentially Expressed Genes (DEGs) based on their promoter sequences using a Convolutional Neural Network (CNN). The approach leverages biological prior knowledge by focusing on promoter regions instead of generic gene features, aiming to capture sequence motifs and regulatory signals that could explain gene expression differences.

---

## Workflow

### 1. Promoter Extraction

- Extract promoter sequences upstream of gene transcription start sites (TSS) from genome annotations (GFF3) and genome sequences (FASTA).
- Typically, a fixed length window upstream (e.g., 1000 bp) is used to capture relevant regulatory regions.
- The script parses the GFF3 file to identify gene coordinates and extracts corresponding promoter sequences from the genome FASTA.

### 2. One-Hot Encoding of Promoter Sequences

- Promoter sequences are converted into one-hot encoded matrices where each nucleotide (A, C, G, T) is represented by a binary vector.
- The output shape is `(number_of_sequences, promoter_length, 4)`, suitable as input for CNN models.
- This encoding preserves spatial information and nucleotide identity, enabling motif detection by convolutional filters.

### 3. CNN Model Training

- A CNN architecture is trained on the one-hot encoded promoter sequences to classify genes as DEG or non-DEG.
- The network typically includes convolutional layers with ReLU activations, max-pooling layers for spatial reduction, dropout for regularization, and fully connected layers ending with a sigmoid output for binary classification.
- Data is split into training and test sets with stratification on class labels.
- Imbalanced classes can be handled by techniques such as SMOTE or class weighting.
- The model is trained with binary cross-entropy loss and evaluated using metrics like ROC AUC, precision, recall, and F1-score.

### 4. Model Prediction and Interpretation

- The trained CNN model can be used to predict DEG status on new promoter sequences.
- Predictions include binary class labels and probabilities.
- For biological interpretability, saliency maps or gradient-based attribution methods highlight important nucleotide positions influencing the model's decision.
- Visualization of these importance scores can guide the discovery of regulatory motifs associated with differential expression.

---

## Requirements

- Python 3.7+
- Biopython
- pandas
- numpy
- scikit-learn
- imblearn (for SMOTE)
- tensorflow (or tensorflow-gpu)
- matplotlib
- seaborn
- gffutils

---

## Directory Structure

├── data/

│ ├── genome.fasta # Reference genome sequences

│ ├── annotation.gff3 # Gene annotation GFF3 file

├── features/

│ ├── gene_features.csv # Optional: gene-level features

│ ├── promoter_sequences.fasta # Extracted promoter sequences

│ ├── promoter_encoded.npy # One-hot encoded promoter sequences

├── models/

│ ├── cnn_promoter_model.h5 # Trained CNN model file

├── scripts/

│ ├── extract_promoters.py # Extract promoters from GFF3+FASTA

│ ├── encode_promoters.py # One-hot encode promoter sequences

│ ├── train_cnn.py # Train CNN classifier on promoters

│ ├── predict_cnn.py # Predict and interpret new sequences

├── README.md


---

## Usage

### Extract Promoters

Run the extraction script with the genome FASTA and GFF3 annotation as input to obtain promoter sequences in FASTA format.

### One-Hot Encoding

Encode the extracted promoter sequences into one-hot format, preparing the input data for the CNN.

### Train CNN

Train the CNN model on the encoded promoters using the DEG labels. Use proper train/test splits and consider data balancing methods.

### Predict New Samples

Use the trained model to predict DEG status for new promoter sequences and generate saliency maps for interpretability.

---

## Notes and Recommendations

- Promoter length choice is critical and may affect model performance. Typical ranges are 500-2000 bp upstream of TSS.
- One-hot encoding preserves nucleotide order and identity, which CNNs exploit for motif learning.
- The CNN architecture and hyperparameters should be tuned based on cross-validation results.
- Saliency or other explainability methods help relate model predictions to biologically meaningful sequence features.
- Integrating other omics data or genomic context may improve model accuracy beyond sequence-only features.

---

## References

- Alipanahi et al., Predicting the sequence specificities of DNA- and RNA-binding proteins by deep learning, Nature Biotechnology, 2015.
- Kelley et al., Basset: learning the regulatory code of the accessible genome with deep convolutional neural networks, Genome Research, 2016.
- Avsec et al., Base-resolution models of transcription-factor binding reveal soft motif syntax, Nature Genetics, 2021.

---

This README outlines the methodology and provides a structured overview of the promoter-based CNN DEG classification pipeline. For detailed usage, see individual scripts in the `scripts/` directory.
