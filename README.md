# ENM-MD多尺度蛋白质动力学模拟框架

## 项目概述
这是一个科研项目，旨在构建弹性网络模型（Elastic Network Model, ENM）和分子动力学（Molecular Dynamics, MD）耦合的多尺度框架，用于解析蛋白质折叠机制。项目实现了核心模块，包括ENM网络构建、多尺度耦合求解器、分析工具等，支持数值算法的准确性和计算效率平衡。

## 项目结构
```
ENM-MD-Coupling/
├── src/                # 源代码
│   ├── enm/            # ENM相关模块
│   ├── coupling/       # 多尺度耦合模块
│   ├── analysis/       # 分析模块
│   ├── io/             # 输入输出模块
│   └── utils/          # 工具模块
├── examples/           # 示例脚本
├── tests/              # 单元测试
├── data/               # 数据文件（如PDB）
├── docs/               # 文档
├── requirements.txt    # 依赖
├── setup.py            # 安装脚本
└── README.md           # 本文档
```

## 技术栈
- Python 3.x
- 核心依赖：numpy, scipy, matplotlib, pandas, scikit-learn, numba, MDAnalysis, biopython, seaborn, pytest
- 优化：Numba JIT加速、SciPy稀疏矩阵、并行计算

## 安装
1. 克隆仓库或下载项目。
2. 安装依赖：
   ```
   pip install -r requirements.txt
   ```
3. （可选）安装项目作为包：
   ```
   python setup.py install
   ```

## 使用示例
运行泛素蛋白（1UBQ）模拟示例：
```
$env:PYTHONPATH = "src"; python examples/1ubq_example.py
```
输出示例：
- B因子相关性验证（与参考数据的Pearson相关系数是否>0.75）
- 前5个B因子值
- 关键残基列表（基于PCA权重）

解释：示例进行ENM模拟、模态分析、B因子计算和关键残基识别。相关性验证可能为False，因为使用简化轨迹；实际应用中需真实MD数据以达到r>0.75。

## 运行测试
验证项目功能和稳定性：
```
$env:PYTHONPATH = "src"; python -m pytest tests/ --verbose
```
测试包括ENM构建、耦合求解、分析验证等，使用1UBQ作为主要测试系统，检查相关性r>0.75和数值稳定性。

## 文档
- `docs/README.md`：概述
- `docs/theory.md`：理论背景和数学公式
- `docs/api_reference.md`：API参考

## 当前状态
项目框架完整，核心功能实现，包括性能优化（如自适应步长、稀疏Hessian）。所有测试通过。进一步工作可包括添加更多示例、集成完整MD模拟或扩展可视化。

如有问题或贡献，请联系作者。
