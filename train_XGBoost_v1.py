import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from sklearn.impute import SimpleImputer
from imblearn.over_sampling import SMOTE
from xgboost import XGBClassifier

# Config
INPUT_FILE = "features/deg_dataset_v2.csv"
TARGET_COL = "DEG"
RANDOM_SEED = 42

# Load dataset
print("Loading dataset...")
df = pd.read_csv(INPUT_FILE)

# Split features and target
X = df.drop(columns=["gene_id", TARGET_COL])
y = df[TARGET_COL]

# Select only numerical features
X = X.select_dtypes(include=["number"])

# Impute missing values
print("Imputing missing values...")
imputer = SimpleImputer(strategy="mean")
X_imputed = pd.DataFrame(imputer.fit_transform(X), columns=X.columns)

# Split into train and test
X_train, X_test, y_train, y_test = train_test_split(
    X_imputed, y, test_size=0.2, stratify=y, random_state=RANDOM_SEED
)

# Apply SMOTE
print("Applying SMOTE...")
smote = SMOTE(random_state=RANDOM_SEED)
X_train_bal, y_train_bal = smote.fit_resample(X_train, y_train)

# Train XGBoost model
print("Training XGBoost model...")
model = XGBClassifier(
    n_estimators=100,
    max_depth=5,
    learning_rate=0.1,
    scale_pos_weight=6.0,  # adjust according to class imbalance
    eval_metric='logloss',
    random_state=RANDOM_SEED
)
model.fit(X_train_bal, y_train_bal)

# Cross-validation
print("Running 5-fold cross-validation...")
cv_scores = cross_val_score(model, X_train_bal, y_train_bal, cv=5, scoring='roc_auc')
print(f"CV ROC AUC scores: {cv_scores}")
print(f"Mean CV ROC AUC: {cv_scores.mean():.3f}")

# Test evaluation
print("Evaluating on test set...")
y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)[:, 1]

print(classification_report(y_test, y_pred))
print("Confusion matrix:")
print(confusion_matrix(y_test, y_pred))
print(f"ROC AUC on test set: {roc_auc_score(y_test, y_prob):.3f}")

# Feature importance plot
importances = pd.Series(model.feature_importances_, index=X.columns).sort_values(ascending=False)

plt.figure(figsize=(8, 5))
sns.barplot(x=importances.values, y=importances.index, palette="magma")
plt.title("Feature Importance - XGBoost")
plt.xlabel("Importance")
plt.tight_layout()
plt.savefig("features/feature_importance_XGBoost.png")
plt.close()

print("Done. Results saved in 'features/' directory.")
