import pytest
import numpy as np
from src.enm.network_builder import ENMNetworkBuilder

def test_adjacency_matrix():
    builder = ENMNetworkBuilder(cutoff=5.0)
    coords = np.array([[0,0,0], [3,0,0], [6,0,0]])
    adj = builder.build_adjacency_matrix(coords)
    assert adj[0,1] == 1.0
    assert adj[0,2] == 0.0
    assert adj[1,2] == 1.0

def test_hessian():
    builder = ENMNetworkBuilder(cutoff=5.0, spring_constant=1.0)
    coords = np.array([[0,0,0], [3,0,0]])
    adj = builder.build_adjacency_matrix(coords)
    H = builder.compute_hessian(coords, adj)
    assert H.shape == (6,6)  # 2 atoms * 3

def test_real_pdb():
    builder = ENMNetworkBuilder()
    univ = builder.load_pdb('data/pdb/1ubq.pdb')
    assert len(univ.atoms) > 0
    coords, adj = builder.generate_topology(univ)
    assert adj.shape[0] == len(coords)
