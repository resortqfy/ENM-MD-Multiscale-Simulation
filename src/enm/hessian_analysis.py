"""模块docstring"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from typing import Tuple

def sparsify_hessian(H: csr_matrix, threshold: float = 1e-5) -> csr_matrix:
    """稀疏化Hessian矩阵
    
    Args:
        H: 输入Hessian矩阵
        threshold: 阈值比例
        
    Returns:
        稀疏化后的矩阵
    """
    max_val = np.max(np.abs(H.data))
    mask = np.abs(H.data) > threshold * max_val
    H.data[~mask] = 0
    H.eliminate_zeros()
    return H

def mode_analysis(H: csr_matrix, num_modes: int = 6) -> Tuple[np.ndarray, np.ndarray]:
    """特征值分解求解低频模式"""
    eigenvalues, eigenvectors = eigsh(H, k=num_modes, which='SM')
    return eigenvalues, eigenvectors

# 示例
# H_sparse = sparsify_hessian(H)
# evals, evecs = mode_analysis(H_sparse)
