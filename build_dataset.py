import pandas as pd

# Paths
GENE_FEATURES = "features/gene_features_v3.csv"
DESEQ_RESULTS = "/home/isagallegor/M.roreri_Side_Quests/data/DEGxML_dataset_base.csv"  
OUTPUT = "features/deg_dataset_v2.csv"

# Parameters to define DEG
PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1

# Load genomic features
print("Reading features...")
features_df = pd.read_csv(GENE_FEATURES)

# Load DESeq2 results and label DEGs
print("Reading DESeq results...")
deseq_df = pd.read_csv(DESEQ_RESULTS)

# Detect gene ID column if needed
if "gene_id" not in deseq_df.columns:
    possible_id_cols = ["Unnamed: 0", "gene", "row.names"]
    for col in possible_id_cols:
        if col in deseq_df.columns:
            deseq_df = deseq_df.rename(columns={col: "gene_id"})
            print(f"Column renamed '{col}' → 'gene_id'")
            break

# Clean gene_id
deseq_df["gene_id"] = deseq_df["gene_id"].str.replace(r"\.m\d+$", "", regex=True)

# Filter necessary columns and create DEG label
required_cols = ["gene_id", "log2FoldChange", "padj"]
if not all(col in deseq_df.columns for col in required_cols):
    raise ValueError(f"File must contain columns: {required_cols}")

deseq_df = deseq_df[required_cols].copy()
deseq_df["DEG"] = ((deseq_df["padj"] < PADJ_THRESHOLD) & 
                   (deseq_df["log2FoldChange"].abs() >= LOG2FC_THRESHOLD)).astype(int)

# Join features and labels
print("Joining datasets...")
merged_df = pd.merge(features_df, deseq_df[["gene_id", "DEG"]], on="gene_id", how="inner")

# Save the final dataset
print(f"Total genes in final dataset: {len(merged_df)}")
merged_df.to_csv(OUTPUT, index=False)
print(f"Complete dataset saved in: {OUTPUT}")
