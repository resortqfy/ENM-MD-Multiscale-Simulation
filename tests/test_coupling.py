import pytest
import numpy as np
from src.coupling.multiscale_solver import MultiscaleSolver
from scipy.sparse import csr_matrix

def test_solver_step():
    mass = np.ones(3)
    gamma = 1.0
    H = csr_matrix(np.eye(3))
    solver = MultiscaleSolver(mass, gamma, H)
    r = np.zeros(3)
    v = np.ones(3)
    dt = 0.1
    F_c = np.zeros(3)
    F_r = np.zeros(3)
    r_new, v_new = solver.step(r, v, dt, F_c, F_r)
    assert np.allclose(r_new, [0.09, 0.09, 0.09], atol=1e-5)
