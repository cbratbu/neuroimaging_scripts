import os
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF

input_path = 'path to data file' # CSV with rows=subjects, columns=lesion burden per ROI
output_path = 'output path' # output directory for W and H matrices

# Data import
tissue_spared = pd.read_csv(input_path, index_col=0)
X = tissue_spared.values.T # Transpose  data so  rows correspond to ROIs as used by NMF

# Run NMF
max_iter = 10000
k = 5

model = NMF(n_components=k, init='nndsvd', max_iter=max_iter, random_state=0)
W = model.fit_transform(X)
H = model.components_

# Calculate reconstruction error
X_approx = np.dot(W, H)
reconstruction_error = np.linalg.norm(X - X_approx, 'fro')
print("Reconstruction error:", reconstruction_error)

# Create W and H outputs
component_tables = {}
W_transposed = W.T
for i in range(k):
    component_table = pd.DataFrame({
        'Region': tissue_spared.columns,
        'Contribution': W_transposed[i, :]
    })
    component_table = component_table.sort_values(by='Contribution', ascending=False).reset_index(drop=True)
    component_tables[f'Component_{i + 1}_values'] = component_table

H_data = pd.DataFrame(H.T, columns=[f'Atom_{i+1}' for i in range(k)])
W_data = pd.DataFrame(W_transposed, columns=tissue_spared.columns)
H_data.to_csv(os.path.join(output_path, "H_data.csv"), index=False)
W_data.to_csv(os.path.join(output_path, "W_data.csv"), index=False)
