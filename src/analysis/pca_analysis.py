import numpy as np
from sklearn.decomposition import PCA

def perform_pca(trajectory: np.ndarray, n_components: int = 10) -> dict:
    """执行PCA分析"""
    pca = PCA(n_components=n_components)
    reduced = pca.fit_transform(trajectory)
    return {'explained_variance': pca.explained_variance_ratio_, 'components': pca.components_}

