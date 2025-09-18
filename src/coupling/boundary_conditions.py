import numpy as np
from typing import List, Union

class BoundaryConditions:
    def __init__(self, fixed_atoms: List[int] = None, periodic: bool = False):
        self.fixed = np.array(fixed_atoms or [])
        self.periodic = periodic

    def apply(self, positions: np.ndarray, forces: np.ndarray, box_size: Union[float, np.ndarray] = None) -> tuple:
        """应用边界条件"""
        if self.periodic and box_size is not None:
            positions %= box_size
        forces[self.fixed] = 0.0
        return positions, forces

# 示例
# bc = BoundaryConditions([0,1,2])
# forces = bc.apply_fixed(forces)
