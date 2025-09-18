import sys
sys.path.append('../src')

from enm.network_builder import ENMNetworkBuilder
from enm.hessian_analysis import sparsify_hessian, mode_analysis

# 类似1ubq的实现
# 假设PDB文件为'spike_ace2.pdb'
builder = ENMNetworkBuilder(cutoff=15.0)  # 调整参数
univ = builder.load_pdb('data/pdb/spike_ace2.pdb')
# ... 继续类似逻辑

