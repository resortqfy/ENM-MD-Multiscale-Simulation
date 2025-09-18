import numpy as np
from sklearn.decomposition import PCA
from typing import Tuple

def enm_pca_comparison(enm_modes: np.ndarray, trajectory: np.ndarray) -> float:
    """计算ENM模式与PCA的相关性"""
    pca = PCA(n_components=enm_modes.shape[1])
    pca.fit(trajectory)
    correlation = np.corrcoef(enm_modes.flatten(), pca.components_.flatten())[0,1]
    return correlation

# 示例
# corr = enm_pca_comparison(evecs, traj)

