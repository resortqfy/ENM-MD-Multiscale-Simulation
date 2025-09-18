import numpy as np
from scipy.sparse import lil_matrix
from typing import Tuple
import MDAnalysis as mda
from numba import jit
from joblib import Parallel, delayed

@jit(nopython=True)
def build_adjacency_matrix(coords: np.ndarray, cutoff: float) -> np.ndarray:
    coords = coords.astype(np.float64)  # 确保float
    n = coords.shape[0]
    adj = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i+1, n):
            dist = np.linalg.norm(coords[i] - coords[j])
            if dist < cutoff:
                adj[i, j] = adj[j, i] = 1.0
    return adj

class ENMNetworkBuilder:
    def __init__(self, cutoff: float = 12.0, spring_constant: float = 1.0):
        self.cutoff = cutoff
        self.k = spring_constant

    def load_pdb(self, pdb_file: str) -> mda.Universe:
        """加载PDB文件并进行预处理"""
        return mda.Universe(pdb_file)

    def generate_topology(self, universe: mda.Universe) -> Tuple[np.ndarray, np.ndarray]:
        """生成弹性网络拓扑"""
        coords = universe.atoms.positions
        adj = build_adjacency_matrix(coords, self.cutoff)
        return coords, adj

    def compute_hessian(self, coords: np.ndarray, adj: np.ndarray) -> lil_matrix:
        """计算Hessian矩阵（稀疏形式）"""
        n = coords.shape[0]
        H = lil_matrix((3*n, 3*n))
        for i in range(n):
            for j in range(n):
                if adj[i, j] > 0 and i != j:
                    diff = coords[i] - coords[j]
                    dist = np.linalg.norm(diff)
                    unit = diff / dist
                    for a in range(3):
                        for b in range(3):
                            val = -self.k * unit[a] * unit[b]
                            H[3*i + a, 3*j + b] = val
                            H[3*j + b, 3*i + a] = val
                            H[3*i + a, 3*i + b] -= val  # 对角元素
                            H[3*j + b, 3*j + a] -= val
        return H

    def build_adjacency_matrix_parallel(self, coords: np.ndarray, n_jobs: int = -1) -> np.ndarray:
        n = coords.shape[0]
        def compute_row(i):
            row = np.zeros(n)
            for j in range(i+1, n):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < self.cutoff:
                    row[j] = 1.0
            return row
        
        rows = Parallel(n_jobs=n_jobs)(delayed(compute_row)(i) for i in range(n))
        adj = np.zeros((n, n))
        for i, row in enumerate(rows):
            adj[i] = row
            adj[:, i] = row  # 对称
        adj += np.eye(n)  # 自连接可选
        return adj

    def build_adjacency_matrix(self, coords: np.ndarray) -> np.ndarray:
        return build_adjacency_matrix(coords, self.cutoff)

# 示例使用
# builder = ENMNetworkBuilder()
# univ = builder.load_pdb('data/pdb/1ubq.pdb')
# coords, adj = builder.generate_topology(univ)
# H = builder.compute_hessian(coords, adj)
