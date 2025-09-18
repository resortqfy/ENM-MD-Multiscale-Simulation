import numpy as np
from sklearn.decomposition import PCA
from typing import List

def bfactor_regression(displacements: np.ndarray, scale_factor: float = 1.0) -> np.ndarray:
    """B因子回归分析"""
    u2 = np.mean(displacements**2, axis=0)
    return (8 * np.pi**2 / 3) * u2 * scale_factor

def pca_weight_analysis(trajectory: np.ndarray, num_components: int = 3) -> List[int]:
    """PCA权重分析识别关键残基"""
    pca = PCA(n_components=num_components)
    pca.fit(trajectory)
    weights = np.abs(pca.components_).sum(axis=0)
    key_residues = np.argsort(weights)[-10:]  # 前10个贡献最大的残基
    return key_residues.tolist()

# 示例
# bfactors = bfactor_regression(displacements)
# keys = pca_weight_analysis(trajectory)

