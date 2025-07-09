import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Dropout, Flatten, Dense, BatchNormalization
from tensorflow.keras.callbacks import EarlyStopping
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
import matplotlib.pyplot as plt

# Config
X_TRAIN_FILE = "promoters/X_train.npy"
X_TEST_FILE = "promoters/X_test.npy"
Y_TRAIN_FILE = "promoters/y_train.npy"
Y_TEST_FILE = "promoters/y_test.npy"
EPOCHS = 30
BATCH_SIZE = 64
RANDOM_SEED = 42

# Load data
print("Loading encoded data...")
X_train = np.load(X_TRAIN_FILE)
X_test = np.load(X_TEST_FILE)
y_train = np.load(Y_TRAIN_FILE)
y_test = np.load(Y_TEST_FILE)

print(f"X_train shape: {X_train.shape}, y_train shape: {y_train.shape}")

# Build CNN model
print("Building model...")
model = Sequential([
    Conv1D(128, 10, activation='relu', input_shape=(X_train.shape[1], 4)),
    BatchNormalization(),
    MaxPooling1D(pool_size=3),
    Dropout(0.3),

    Conv1D(64, 5, activation='relu'),
    MaxPooling1D(pool_size=2),
    Dropout(0.3),

    Flatten(),
    Dense(64, activation='relu'),
    Dropout(0.4),
    Dense(1, activation='sigmoid')
])

model.compile(optimizer='adam',
              loss='binary_crossentropy',
              metrics=['accuracy', tf.keras.metrics.AUC(name='auc')])

model.summary()

# Training
print("Training model...")
early_stop = EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True)

history = model.fit(
    X_train, y_train,
    validation_split=0.2,
    epochs=EPOCHS,
    batch_size=BATCH_SIZE,
    callbacks=[early_stop],
    verbose=1
)

# Evaluation
print("Evaluating on test set...")
y_pred_prob = model.predict(X_test).flatten()
y_pred = (y_pred_prob > 0.5).astype(int)

print("Classification report:")
print(classification_report(y_test, y_pred))
print("Confusion matrix:")
print(confusion_matrix(y_test, y_pred))
print(f"ROC AUC: {roc_auc_score(y_test, y_pred_prob):.3f}")

# Save model
model.save("promoters/cnn_promoter_model.h5")

# Plot training
plt.figure(figsize=(10, 4))
plt.plot(history.history['loss'], label='Train loss')
plt.plot(history.history['val_loss'], label='Val loss')
plt.legend()
plt.title("Training and validation loss")
plt.tight_layout()
plt.savefig("promoters/cnn_training_loss.png")
plt.close()
