import numpy as np
from scipy.sparse import csr_matrix
from typing import Callable, Tuple
from .adaptive_integrator import AdaptiveIntegrator

class MultiscaleSolver:
    def __init__(self, mass: np.ndarray, gamma: float, H: csr_matrix, adaptive: bool = False):
        self.M = mass
        self.gamma = gamma
        self.H = H
        self.adaptive = adaptive
        self.integrator = AdaptiveIntegrator() if adaptive else None

    def step(self, r: np.ndarray, v: np.ndarray, dt: float, F_coupling: np.ndarray, F_random: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """单步积分耦合动力学方程"""
        acc = (-self.gamma * v - self.H @ r + F_coupling + F_random) / self.M
        v_new = v + acc * dt
        r_new = r + v_new * dt
        return r_new, v_new

    def simulate(self, r0: np.ndarray, v0: np.ndarray, t_max: float, dt: float, force_fn: Callable) -> np.ndarray:
        """运行模拟"""
        r, v = r0.copy(), v0.copy()
        trajectory = [r.copy()]
        t = 0
        while t < t_max:
            F_coup, F_rand = force_fn(r, v, t)
            acc = (-self.gamma * v - self.H @ r + F_coup + F_rand) / self.M
            if self.adaptive:
                dt = self.integrator.get_step(acc, v)
            r, v = self.step(r, v, dt, F_coup, F_rand)
            trajectory.append(r.copy())
            t += dt
        return np.array(trajectory)
