# 基于弹性网络模型的蛋白质关键残基识别与动力学多尺度模拟
# Identification of Key Protein Residues and Multiscale Dynamic Simulation Based on an Elastic Network Model

**项目负责人 (PI)**: 乔斐远 (Feiyuan Qiao)
**团队成员 (Team)**: 柳卓炜 (Zhuowei Liu), 王嘉悦 (Jiayue Wang)
**指导教师 (Supervisor)**: 蔡勇勇 教授 (Prof. Yongyong Cai)
**单位**: 北京师范大学 数学科学学院 (Beijing Normal University, School of Mathematical Sciences)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
> **注意**: 上面的 DOI 徽章是占位符。当你项目完成，按照我之前提到的方法在 Zenodo 上发布后，替换成你自己的 DOI。

---

## 摘要 (Abstract)

本项目旨在通过计算机模拟方法，以微分方程为核心，探索构建弹性网络模型（ENM）与分子动力学（MD）的多尺度耦合框架。我们结合分子动力学领域的现有研究成果，针对特定蛋白质分子（如泛素蛋白 1UBQ）开展模拟，以验证模型的可靠性和有效性，为蛋白质折叠机制解析提供新的数学模型。

This project aims to construct a multiscale coupling framework combining the Elastic Network Model (ENM) and Molecular Dynamics (MD), driven by differential equations. By simulating specific proteins like Ubiquitin (1UBQ), we validate the model's reliability and effectiveness, offering a new mathematical paradigm for analyzing protein folding mechanisms.

## 项目结构 (Project Structure)

```
.
├── data/               # 存放原始 PDB 文件和实验数据
├── docs/               # 存放项目报告和最终论文
├── results/            # 存放模拟结果、图表和分析数据
├── scripts/            # 存放主要的 Python 模拟和分析脚本
│   ├── 01_preprocess.py
│   ├── 02_run_simulation.py
│   ├── 03_analyze_rmsf.py
│   └── ...
├── .gitignore          # Git 忽略文件配置
├── LICENSE             # MIT 许可证
└── README.md           # 本说明文件
```

## 环境要求 (Requirements)

*   Python 3.8+
*   GROMACS
*   NumPy
*   SciPy
*   (其他你可能用到的库)

为了方便复现，请使用以下命令安装所需 Python 依赖：
```bash
pip install -r requirements.txt
```

## 使用方法 (Usage)

1.  **准备工作**: 将蛋白质的 PDB 文件（如 `1ubq.pdb`）放入 `data/` 文件夹。

2.  **预处理与模拟**:
    ```bash
    # 运行 GROMACS 预处理和全原子 MD 模拟 (详细步骤请参考 docs/)
    # ...
    ```

3.  **运行多尺度模型与分析**:
    ```bash
    # 运行 ENM-MD 耦合模拟
    python scripts/02_run_simulation.py --input data/1ubq.pdb --output results/trajectory.xtc

    # 分析 RMSF 并识别关键残基
    python scripts/03_analyze_rmsf.py --input results/trajectory.xtc --output results/key_residues.csv
    ```

## 如何引用 (Citation)

如果我们的工作对您有帮助，请引用我们的论文和这份代码。
If you find our work useful in your research, please consider citing our paper and this repository.

```bibtex
@article{qiao2026protein,
  title   = {Identification of Key Protein Residues and Multiscale Dynamic Simulation Based on an Elastic Network Model},
  author  = {Qiao, Feiyuan and Liu, Zhuowei and Wang, Jiayue and Cai, Yongyong},
  journal = {Journal of Computational Biology},
  year    = {2026},
  url     = {https://github.com/[你的用户名]/[你的仓库名]}
}
```
> **注意**: 上面的 BibTeX 是一个模板，请在论文发表后更新为实际信息。

## 致谢 (Acknowledgements)

本研究受北京师范大学本科生科研训练与创新创业项目资助。
This work was supported by the Undergraduate Research Training Program and Innovation and Entrepreneurship Project of Beijing Normal University.
