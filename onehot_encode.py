import pandas as pd
import numpy as np
import os
from sklearn.model_selection import train_test_split
import pickle

# Config
CSV_IN = "promoters/promoters.csv"
X_OUT = "promoters/X_onehot.npy"
Y_OUT = "promoters/y_labels.npy"
SPLIT = True
X_TRAIN_OUT = "promoters/X_train.npy"
X_TEST_OUT = "promoters/X_test.npy"
Y_TRAIN_OUT = "promoters/y_train.npy"
Y_TEST_OUT = "promoters/y_test.npy"
RANDOM_SEED = 42

# One-hot encoding 
def onehot_encode_seq(seq):
    mapping = {"A": [1, 0, 0, 0],
               "C": [0, 1, 0, 0],
               "G": [0, 0, 1, 0],
               "T": [0, 0, 0, 1],
               "N": [0, 0, 0, 0]}  # N represents ambiguous nucleotides 

    seq = seq.upper()
    return np.array([mapping.get(nt, [0, 0, 0, 0]) for nt in seq])

# Load promoters file
print("Loading promoter sequences...")
df = pd.read_csv(CSV_IN)
sequences = df["promoter_seq"]
labels = df["DEG"].values.astype(int)

# Encode every sequence
print("Encoding sequences...")
onehot_seqs = np.array([onehot_encode_seq(seq) for seq in sequences])

# Validate dimentions
print(f"Encoded shape: {onehot_seqs.shape}  (num_sequences, seq_len, 4)")

# Save and split
if SPLIT:
    print("Splitting train/test sets...")
    X_train, X_test, y_train, y_test = train_test_split(
        onehot_seqs, labels, test_size=0.2, random_state=RANDOM_SEED, stratify=labels
    )
    np.save(X_TRAIN_OUT, X_train)
    np.save(X_TEST_OUT, X_test)
    np.save(Y_TRAIN_OUT, y_train)
    np.save(Y_TEST_OUT, y_test)
    print("Saved train/test splits.")
else:
    np.save(X_OUT, onehot_seqs)
    np.save(Y_OUT, labels)
    print("Saved all data as X and y.")
