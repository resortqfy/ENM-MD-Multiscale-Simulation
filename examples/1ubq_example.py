import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


from enm.network_builder import ENMNetworkBuilder
from enm.hessian_analysis import sparsify_hessian, mode_analysis
from analysis.residue_analysis import bfactor_regression, pca_weight_analysis
from enm.elastic_dynamics import ENMDynamics
from utils.constants import GAMMA
from scipy.sparse import csr_matrix
import pandas as pd
from analysis.validation import validate_bfactor_correlation
import numpy as np
import MDAnalysis as mda

# 加载PDB
builder = ENMNetworkBuilder()
univ = builder.load_pdb('data/pdb/1ubq.pdb')
coords, adj = builder.generate_topology(univ)
H = builder.compute_hessian(coords, adj)

# Hessian to csr
H_sparse = sparsify_hessian(H.tocsr())

# 模态
evals, evecs = mode_analysis(H_sparse)

# 动力学
masses = np.ones(len(coords) * 3)
dynamics = ENMDynamics(H_sparse, masses)
r0 = coords.flatten()
v0 = np.zeros_like(r0)
traj = dynamics.solve(r0, v0, 0.01, 100)

# 提取参考B因子 (per atom, then average per residue)
atoms = univ.select_atoms('protein')
ref_bfactors_atoms = atoms.tempfactors  # B因素
residues = np.unique(atoms.resids)
ref_bfactors = np.array([np.mean(ref_bfactors_atoms[atoms.resids == res]) for res in residues])

# 调整计算bfactors为per-residue
n_res = len(residues)
n_frames = traj.shape[0]
traj_res = np.zeros((n_frames, n_res, 3))
for i, res in enumerate(residues):
    res_atoms = univ.select_atoms(f'resid {res}')
    traj_res[:, i] = traj.reshape(n_frames, -1, 3)[:, res_atoms.ix].mean(axis=1)

displacements = traj_res - traj_res.mean(axis=0)
bfactors = bfactor_regression(displacements.reshape(n_frames, n_res * 3).T.reshape(n_res, 3 * n_frames).T)  # 调整形状

is_valid = validate_bfactor_correlation(bfactors, ref_bfactors)
print("B因子相关性验证:", is_valid)

keys = pca_weight_analysis(traj.reshape(-1, len(coords), 3).mean(axis=2))

print("B因素:", bfactors[:5])
print("关键残基:", keys)
