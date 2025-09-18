import numpy as np
from numba import jit

@jit(nopython=True)
def compute_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """计算RMSD"""
    return np.sqrt(np.mean((a - b)**2))

