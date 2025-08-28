import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from scipy.linalg import svd

# set paths
tissue_data_path = 'path to lesion burden data' # CSV with rows of patients, columns of ROIs (IDs in first column)

# Data imports
tissue_spared = pd.read_csv(tissue_data_path, index_col=0)
X = tissue_spared.values.T

# random matrix with the same mean and variance as the original data matrix
mean_X = np.mean(X)
std_X = np.std(X)
random_matrix = np.random.normal(mean_X, std_X, X.shape)

# Ensure non-negativity
random_matrix[random_matrix < 0] = 0

# Range of k values to try
k_values = range(1, 21)
nmf_reconstruction_errors = []
svd_reconstruction_errors = []
svd_random_reconstruction_errors = []

# NMF
for k in k_values:
    model = NMF(n_components=k, init='nndsvd', max_iter=10000, random_state=0)
    W = model.fit_transform(X)
    H = model.components_
    X_approx = np.dot(W, H)
    reconstruction_error = np.linalg.norm(X - X_approx, 'fro')
    nmf_reconstruction_errors.append(reconstruction_error)
    print(f'NMF: k={k}, reconstruction error={reconstruction_error}')

# SVD
for k in k_values:
    U, s, Vt = svd(X, full_matrices=False)
    Uk = U[:, :k]
    Sk = np.diag(s[:k])
    Vk = Vt[:k, :]
    X_approx_svd = np.dot(Uk, np.dot(Sk, Vk))
    reconstruction_error = np.linalg.norm(X - X_approx_svd, 'fro')
    svd_reconstruction_errors.append(reconstruction_error)
    print(f'SVD: k={k}, reconstruction error={reconstruction_error}')

# SVD on random matrix
for k in k_values:
    U, s, Vt = svd(random_matrix, full_matrices=False)
    Uk = U[:, :k]
    Sk = np.diag(s[:k])
    Vk = Vt[:k, :]
    X_approx_svd_random = np.dot(Uk, np.dot(Sk, Vk))
    reconstruction_error = np.linalg.norm(random_matrix - X_approx_svd_random, 'fro')
    svd_random_reconstruction_errors.append(reconstruction_error)
    print(f'SVD (Random Matrix): k={k}, reconstruction error={reconstruction_error}')

# Plot the reconstruction error for NMF, SVD, and SVD on random data against k values
plt.figure(figsize=(10, 6))
plt.plot(k_values, nmf_reconstruction_errors, marker='o', linestyle='--', label='NMF')
plt.plot(k_values, svd_reconstruction_errors, marker='x', linestyle='-', label='SVD')
plt.plot(k_values, svd_random_reconstruction_errors, marker='s', linestyle=':', label='SVD (Random Matrix)')

plt.title('Reconstruction Error vs. Number of Components (k)')
plt.xlabel('Number of Components (k)')
plt.ylabel('Reconstruction Error (RMS)')
plt.xticks(ticks=k_values)
plt.legend()
plt.grid(True)
plt.show()
