import pytest
from src.enm.network_builder import ENMNetworkBuilder
from src.analysis.validation import validate_correlation
from src.enm.elastic_dynamics import ENMDynamics
from src.utils.math_utils import compute_rmsd
import numpy as np

def test_integration():
    # 简化集成测试
    builder = ENMNetworkBuilder()
    # 假设数据
    computed = np.random.randn(10)
    ref = computed * 0.8 + np.random.randn(10)*0.1
    assert validate_correlation(computed, ref)

def test_full_simulation():
    # 假设数据
    H = np.eye(3)
    masses = np.ones(3)
    dynamics = ENMDynamics(H, masses)
    r0 = np.zeros(3)
    v0 = np.ones(3)
    traj = dynamics.solve(r0, v0, 0.1, 10)
    assert traj.shape[0] == 11
    rmsd = compute_rmsd(traj[0], traj[-1])
    assert rmsd > 0
