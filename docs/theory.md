# 理论背景

## ENM模型
弹性网络模型（ENM）将蛋白质简化为连接残基的弹性网络。势能函数：
```
V_ENM = 0.5 * k * Σ (r_ij - r_ij^0)^2  (对所有连接对)
```
其中k为弹性常数，r_ij为当前距离，r_ij^0为平衡距离。Hessian矩阵从此导出，用于模态分析。

## 多尺度耦合
ENM-MD耦合通过方程：
```
M * d²r/dt² + γ * dr/dt + H * r = F_coupling + F_random
```
整合粗粒ENM和细粒MD。自适应步长基于曲率：
```
kappa = ||d²r/dt²|| / (||dr/dt|| + epsilon)
dt = dt_max / (1 + alpha * kappa)
```

## 分析方法
- B因子回归：B = (8π²/3) * <u²> * scale
- PCA：识别主导模式贡献残基
- 验证：Pearson相关系数r > 0.75与实验数据

参考文献：
- Vijay-Kumar et al. (1987) on Ubiquitin structure
- ENM理论：Tirion (1996)
