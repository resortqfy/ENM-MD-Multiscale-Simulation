import MDAnalysis as mda
import MDAnalysis.analysis.rms
import numpy as np
import argparse
import os

def calculate_per_residue_rmsd(tpr_file, xtc_file, start_res, end_res, fit_select, calc_type, output_prefix="rmsd_res"):
    """
    计算蛋白质中指定残基范围的逐残基 RMSD。

    Args:
        tpr_file (str): GROMACS .tpr 文件路径。
        xtc_file (str): 轨迹文件路径 (.xtc, .trr, .dcd 等)。
        start_res (int): 起始残基编号（包含）。
        end_res (int): 结束残基编号（包含）。
        fit_select (str): 用于最小二乘拟合（对齐）的原子选择字符串（例如："protein and backbone"）。
        calc_type (str): 计算 RMSD 的原子类型 ('bb' 表示主链, 'whole' 表示所有原子)。
        output_prefix (str): 输出文件名的前缀。
    """

    print(f"--- MDAnalysis 逐残基 RMSD 计算脚本 ---")
    print(f"加载轨迹：{xtc_file} (拓扑：{tpr_file})")
    print(f"残基范围：{start_res} - {end_res}")
    print(f"对齐组：'{fit_select}'")
    print(f"计算组：每个残基的 '{'主链' if calc_type == 'bb' else '所有原子'}'")
    print("-------------------------------------")

    # 检查文件是否存在
    if not os.path.exists(tpr_file):
        print(f"错误：未找到 TPR 文件 '{tpr_file}'。")
        return
    if not os.path.exists(xtc_file):
        print(f"错误：未找到轨迹文件 '{xtc_file}'.")
        return

    try:
        # 加载宇宙（模拟系统）
        u = mda.Universe(tpr_file, xtc_file)
        print(f"成功加载宇宙，包含 {u.select_atoms('all').n_atoms} 个原子， {u.trajectory.n_frames} 帧。")

        # 获取参考结构 (通常是第一帧)
        reference = u.select_atoms("all")
        u.trajectory[0] # 确保宇宙指针在第一帧

        # 构建需要计算 RMSD 的各个残基的原子选择列表
        residue_selections = []
        residue_numbers = []
        print(f"构建残基选择字符串...")
        for res_i in range(start_res, end_res + 1):
            if calc_type == 'bb':
                selection_string = f"resid {res_i} and backbone"
            else: # whole residue
                selection_string = f"resid {res_i}"

            # 检查选择是否有效且非空
            selected_atoms = u.select_atoms(selection_string)
            if selected_atoms.n_atoms == 0:
                print(f"警告：残基 {res_i} 的选择 '{selection_string}' 未选中任何原子。跳过此残基。")
                continue # 跳过当前残基

            residue_selections.append(selection_string)
            residue_numbers.append(res_i) # 记录对应残基的编号

        if not residue_selections:
            print("错误：没有找到任何有效的残基进行计算。请检查残基范围和类型选择是否正确。")
            return

        print(f"将对 {len(residue_selections)} 个有效残基 ({residue_numbers[0]}-{residue_numbers[-1]}) 进行计算。")

        # 创建 RMSD 分析对象
        # select 参数用于指定对齐组
        # groupselections 参数用于指定计算 RMSD 的组列表
        R = MDAnalysis.analysis.rms.RMSD(
            u,
            reference=reference, # 参考结构
            select=fit_select,   # 对齐组
            groupselections=residue_selections # 计算组列表
        )

        print("开始计算 RMSD...")
        # 运行分析
        R.run()

        print("RMSD 计算完成。正在保存结果...")

        # 提取结果
        # R.results.rmsd 是一个 NumPy 数组
        # 第一列是时间 (皮秒), 之后每一列对应 groupselections 中的一个选择
        results = R.results.rmsd

        # 保存每个残基的 RMSD 到单独的文件
        # results[:, 0] 是时间
        # results[:, j] 是 groupselections[j-1] 对应的 RMSD (索引从1开始对应列)
        for i, res_select_string in enumerate(residue_selections):
            res_num = residue_numbers[i] # 获取对应的残基编号

            # 提取当前残基的 RMSD 数据 (对应 results 数组的 i+1 列)
            current_residue_rmsd = results[:, i + 1]

            # 将时间从皮秒转换为纳秒
            time_ns = results[:, 0] / 1000.0

            # 构建输出文件名
            output_file = f"{output_prefix}_res{res_num}_{calc_type}.xvg"

            # 准备写入的数据 (时间, RMSD)
            data_to_save = np.column_stack((time_ns, current_residue_rmsd))

            # 写入 .xvg 文件
            # 添加一些 GROMACS .xvg 格式的注释
            header = f"""# Per-residue RMSD for residue {res_num} ({calc_type} atoms)
# Aligned to '{fit_select}'
@ title "RMSD of Residue {res_num} ({calc_type})"
@ xaxis label "Time (ns)"
@ yaxis label "RMSD ($\\AA$)"
@ type xy
"""
            np.savetxt(output_file, data_to_save, fmt='%.6f', header=header, comments='')
            print(f"已保存残基 {res_num} 的 RMSD 到 '{output_file}'")

        print("\n脚本执行完毕。")

    except Exception as e:
        print(f"\n发生错误：{e}")
        import traceback
        traceback.print_exc()


# --- 主程序入口 ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="使用 MDAnalysis 计算蛋白质逐残基 RMSD。")
    parser.add_argument("-s", "--tpr", required=True, help="GROMACS .tpr 文件路径。")
    parser.add_argument("-f", "--xtc", required=True, help="轨迹文件路径 (.xtc, .trr, .dcd 等)。")
    parser.add_argument("--start-res", type=int, required=True, help="起始残基编号（包含）。")
    parser.add_argument("--end-res", type=int, required=True, help="结束残基编号（包含）。")
    parser.add_argument("--fit-select", required=True, help="用于最小二乘拟合（对齐）的原子选择字符串（例如：'protein and backbone'）。")
    parser.add_argument("--calc-type", choices=['bb', 'whole'], required=True, help="计算 RMSD 的原子类型 ('bb' 表示主链, 'whole' 表示所有原子)。")
    parser.add_argument("--output-prefix", default="rmsd_res", help="输出文件名的前缀（默认为 'rmsd_res'）。")

    args = parser.parse_args()

    calculate_per_residue_rmsd(
        args.tpr,
        args.xtc,
        args.start_res,
        args.end_res,
        args.fit_select,
        args.calc_type,
        args.output_prefix
    )

