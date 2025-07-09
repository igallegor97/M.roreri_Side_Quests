import numpy as np
from tensorflow.keras.models import load_model
from Bio import SeqIO
from scripts.encode_promoters import one_hot_encode_seq  
import matplotlib.pyplot as plt

# Config
FASTA_FILE = "promoters/new_promoters.fasta"   # File with simple samples
MODEL_FILE = "promoters/cnn_promoter_model.h5"
PROMOTER_LENGTH = 1000  # Must match the one used in training
OUTPUT_PRED = "promoters/new_predictions.csv"

# Load model
print("Loading trained model...")
model = load_model(MODEL_FILE)

# Read and encode sequences
print("Reading and encoding new sequences...")
records = list(SeqIO.parse(FASTA_FILE, "fasta"))
X_new = np.array([one_hot_encode_seq(str(rec.seq).upper(), PROMOTER_LENGTH) for rec in records])

# Predict
print("Predicting...")
y_prob = model.predict(X_new).flatten()
y_pred = (y_prob > 0.5).astype(int)

# Save results
import pandas as pd
results_df = pd.DataFrame({
    "gene_id": [rec.id for rec in records],
    "DEG_prediction": y_pred,
    "probability": y_prob
})
results_df.to_csv(OUTPUT_PRED, index=False)
print(f"Predictions saved to {OUTPUT_PRED}")
