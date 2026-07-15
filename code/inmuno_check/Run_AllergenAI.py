import sys,os
import numpy as np
import pandas as pd
import keras
import tensorflow as tf

input_file = sys.argv[1] # e.g., Cupin.txt
idx_file = sys.argv[2]   # e.g., Cupin_idx.txt
output_file = sys.argv[3] # raw report output

# Load model (make sure finalmodel.h5 is in the same directory or adjust path)
model_path = os.environ.get("ALLERGENAI_MODEL_PATH", "./finalmodel.h5")
model = keras.models.load_model(model_path)

# Load the input file directly into a numpy array (much more memory efficient than pandas)
# The file has N*1000 rows and 20 columns of space-separated 0s and 1s
arr = np.loadtxt(input_file, dtype=np.float32)
protein = arr.reshape(-1, 1000, 20)

pred = model.predict(protein)

# Map back to sequence names using idx_file
# idx_file contains: 1 \t num \n
idx_df = pd.read_csv(idx_file, sep='\t', header=None, names=["c", "num"])

# Prepare output dataframe
results = []
for i in range(len(pred)):
    seq_name = idx_df.iloc[i]["num"]
    p_non_allergen = pred[i][0]
    p_allergen = pred[i][1]
    
    # Check if prediction output matches expected shape
    results.append({
        "Sequence": seq_name,
        "P_non_allergen": p_non_allergen,
        "P_allergen": p_allergen,
        "Prediction": "Allergen" if p_allergen > p_non_allergen else "Non-allergen"
    })

df = pd.DataFrame(results)
df.to_csv(output_file, index=False)
