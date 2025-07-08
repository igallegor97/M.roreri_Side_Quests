import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from sklearn.impute import SimpleImputer
from imblearn.over_sampling import SMOTE
import shap

# CONFIG
INPUT_FILE = "features/deg_dataset_v2.csv"
TARGET_COL = "DEG"
RANDOM_SEED = 42

# Load Dataset
print("Loading dataset...")
df = pd.read_csv(INPUT_FILE)

# Selecting predictrs (x) and labels (y)
X = df.drop(columns=["gene_id", TARGET_COL])
y = df[TARGET_COL]

# Only numeric columns
X = X.select_dtypes(include=["number"])

# Split train/test
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=RANDOM_SEED
)

# Verify columns with NaNs
print("NaNs in X_train before imputation:")
print(X_train.isna().sum())

# Delete columns with NaNs
X_train = X_train.dropna(axis=1, how="all")
X_test = X_test[X_train.columns]  

# Impute missing values with mean
imputer = SimpleImputer(strategy="mean")
X_train = pd.DataFrame(imputer.fit_transform(X_train), columns=X_train.columns)
X_test = pd.DataFrame(imputer.transform(X_test), columns=X_train.columns)

# Verify resulting NaNs 
assert X_train.isna().sum().sum() == 0, "NaNs in X_train"
assert X_test.isna().sum().sum() == 0, "NaNs in X_test"

# SMOTE class balance
print("Applying SMOTE to balance training data...")
smote = SMOTE(random_state=RANDOM_SEED)
X_train_bal, y_train_bal = smote.fit_resample(X_train, y_train)

# Model training
print("Training RF model...")
model = RandomForestClassifier(
    n_estimators=100,
    random_state=RANDOM_SEED,
    class_weight="balanced"
)
model.fit(X_train_bal, y_train_bal)

# Eval
print("Evaluating model in test set...")
y_pred = model.predict(X_test)
print(classification_report(y_test, y_pred))

print("Confusion matrix:")
print(confusion_matrix(y_test, y_pred))

# AUC
y_prob = model.predict_proba(X_test)[:, 1]
auc_score = roc_auc_score(y_test, y_prob)
print(f"ROC AUC: {auc_score:.3f}")

# Feature importance
print("Calculating feature importances...")
importances = pd.Series(
    model.feature_importances_,
    index=X_train.columns
).sort_values(ascending=False)

print(importances)

# Save plot
plt.figure(figsize=(8, 5))
sns.barplot(x=importances.values, y=importances.index, palette="viridis")
plt.title("Feature importance - Random Forest")
plt.xlabel("Importance")
plt.tight_layout()
plt.savefig("features/feature_importance_RF.png")
plt.close()

# SHAP interpretability
print("Calculating SHAP values...")
explainer = shap.Explainer(model, X_train)
shap_values = explainer(X_test)

shap.summary_plot(shap_values, X_test, plot_type="bar", show=False)
plt.tight_layout()
plt.savefig("features/shap_summary_RF.png")
plt.close()

print("Done. Results saved in 'features/' directory.")
