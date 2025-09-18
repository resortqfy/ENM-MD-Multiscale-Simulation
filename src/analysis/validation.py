import numpy as np
from scipy.stats import pearsonr

def validate_correlation(computed: np.ndarray, reference: np.ndarray, threshold: float = 0.75) -> bool:
    """验证相关性"""
    corr = np.corrcoef(computed, reference)[0,1]
    return corr > threshold

def check_stability(energy: np.ndarray) -> bool:
    """检查能量守恒"""
    return np.std(energy) < 1e-3  # 简化检查

def validate_bfactor_correlation(computed: np.ndarray, experimental: np.ndarray, threshold: float = 0.75) -> bool:
    """验证B因子相关性"""
    corr, _ = pearsonr(computed, experimental)
    return corr > threshold
