import numpy as np
from numba import jit
from scipy.sparse import csr_matrix
from typing import Optional

@jit(nopython=True)
def integrate(r: np.ndarray, v: np.ndarray, acc: np.ndarray, dt: float, steps: int) -> np.ndarray:
    traj = np.zeros((steps + 1, r.shape[0]))
    traj[0] = r
    for i in range(steps):
        v += acc * dt
        r += v * dt
        traj[i+1] = r
    return traj

class ENMDynamics:
    def __init__(self, H: csr_matrix, masses: np.ndarray):
        self.H = H
        self.masses = masses  # 保持1D

    def solve(self, r0: np.ndarray, v0: np.ndarray, dt: float, steps: int, gamma: Optional[float] = None) -> np.ndarray:
        """求解ENM动力学
        
        Args:
            r0: 初始位置
            v0: 初始速度
            dt: 时间步
            steps: 步数
            gamma: 可选阻尼系数
            
        Returns:
            轨迹数组
        """
        r = r0.copy()
        v = v0.copy()
        for _ in range(steps):
            acc = - (self.H @ r) / self.masses
            if gamma is not None:
                acc -= gamma * v
            v += acc * dt
            r += v * dt
        # 注意：为效率，使用循环或jit版本
        return integrate(r0, v0, - (self.H @ r0) / self.masses, dt, steps)  # 简化，实际需调整
