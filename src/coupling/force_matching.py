import numpy as np
from sklearn.linear_model import LinearRegression

def force_matching(md_forces: np.ndarray, enm_displacements: np.ndarray) -> np.ndarray:
    """力匹配算法"""
    model = LinearRegression()
    model.fit(enm_displacements, md_forces)
    matched_forces = model.predict(enm_displacements)
    return matched_forces

