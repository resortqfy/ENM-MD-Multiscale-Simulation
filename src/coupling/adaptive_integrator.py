import numpy as np
from numba import jit
from typing import Any

@jit(nopython=True)
def compute_curvature(acc: np.ndarray, vel: np.ndarray, epsilon: float = 1e-6) -> float:
    """计算局部曲率"""
    return np.linalg.norm(acc) / (np.linalg.norm(vel) + epsilon)

@jit(nopython=True)
def adaptive_step(dt_max: float, kappa: float, alpha: float = 1.0) -> float:
    """动态调整步长"""
    return dt_max / (1 + alpha * kappa)

class AdaptiveIntegrator:
    def __init__(self, dt_max: float = 0.01, alpha: float = 1.0):
        self.dt_max = dt_max
        self.alpha = alpha

    def get_step(self, acc: np.ndarray, vel: np.ndarray) -> float:
        kappa = compute_curvature(acc, vel)
        return adaptive_step(self.dt_max, kappa, self.alpha)
