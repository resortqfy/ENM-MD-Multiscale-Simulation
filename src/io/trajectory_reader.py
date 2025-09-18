import MDAnalysis as mda
import numpy as np

def read_trajectory(traj_file: str, top_file: str) -> np.ndarray:
    """读取模拟轨迹"""
    u = mda.Universe(top_file, traj_file)
    coords = np.array([ts.positions for ts in u.trajectory])
    return coords

