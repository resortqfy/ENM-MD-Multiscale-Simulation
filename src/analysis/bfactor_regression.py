import numpy as np
from scipy.stats import linregress

def bfactor_correlation(computed_b: np.ndarray, exp_b: np.ndarray) -> float:
    """计算B因子相关性"""
    slope, intercept, r_value, _, _ = linregress(computed_b, exp_b)
    return r_value

# 示例
# r = bfactor_correlation(computed, experimental)

