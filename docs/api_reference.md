# API参考

## enm模块
- ENMNetworkBuilder：构建ENM网络
  - load_pdb(pdb_file)：加载PDB
  - generate_topology(universe)：生成坐标和邻接矩阵
  - compute_hessian(coords, adj)：计算稀疏Hessian

- hessian_analysis
  - sparsify_hessian(H, threshold)：稀疏化矩阵
  - mode_analysis(H, num_modes)：特征分解

- elastic_dynamics.ENMDynamics：ENM动力学求解
  - solve(r0, v0, dt, steps, gamma)：运行模拟

## coupling模块
- multiscale_solver.MultiscaleSolver：多尺度求解
  - step(r, v, dt, F_coupling, F_random)：单步积分
  - simulate(r0, v0, t_max, dt, force_fn)：完整模拟

- adaptive_integrator.AdaptiveIntegrator：自适应步长
  - get_step(acc, vel)：计算动态dt

- boundary_conditions.BoundaryConditions：边界处理
  - apply(positions, forces, box_size)：应用固定/周期边界

- force_matching.force_matching(md_forces, enm_displacements)：力匹配

## analysis模块
- residue_analysis
  - bfactor_regression(displacements, scale_factor)：计算B因子
  - pca_weight_analysis(trajectory, num_components)：识别关键残基

- bfactor_regression.bfactor_correlation(computed_b, exp_b)：相关性计算

- pca_analysis.perform_pca(trajectory, n_components)：执行PCA

- validation
  - validate_correlation(computed, reference, threshold)：验证相关性
  - validate_bfactor_correlation(computed, experimental, threshold)：B因子验证
  - check_stability(energy)：检查能量守恒

## io模块
- pdb_handler
  - parse_pdb(pdb_file)：解析PDB
  - mda_load(pdb_file)：加载MDAnalysis Universe

- trajectory_reader.read_trajectory(traj_file, top_file)：读取轨迹

- output_writer
  - write_pdb(coords, output_file, template)：写入PDB轨迹
  - write_csv(data, output_file)：写入CSV

## utils模块
- math_utils.compute_rmsd(a, b)：计算RMSD

- plotting
  - plot_bfactors(bfactors, output_file)：绘制B因子
  - visualize_modes(coords, modes, index)：可视化模态（需py3Dmol）

- constants：物理常数如BOLTZMANN, SPRING_CONSTANT等
