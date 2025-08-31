# 文件名: visualize_stiffness.py

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors # 新增导入，用于创建自定义色图
import seaborn as sns

def _create_custom_colormap(vmin, vmax):
    """
    创建一个自定义的蓝-黄-红色图。
    - 负值: 蓝色
    - 零值: 黄色
    - 正值: 红色
    """
    # 如果所有值都为正，则创建从黄到红的色图
    if vmin >= 0:
        print("Data is all positive. Creating a yellow-to-red colormap.")
        return mcolors.LinearSegmentedColormap.from_list(
            "custom_yellow_red", ["yellow", "#FF3333", "red"]
        )
    
    # 如果所有值都为负，则创建从蓝到黄的色图
    if vmax <= 0:
        print("Data is all negative. Creating a blue-to-yellow colormap.")
        return mcolors.LinearSegmentedColormap.from_list(
            "custom_blue_yellow", ["blue", "#FFFF99", "yellow"]
        )

    # --- 核心逻辑：创建以0为中心的分段色图 ---
    # 数据横跨正负值。我们需要计算 '0' 在整个数据范围 [vmin, vmax] 中所处的位置比例。
    # 例如，如果 vmin=-10, vmax=30，那么0点的位置就在 (-(-10))/(30 - (-10)) = 10/40 = 0.25
    zero_norm = -vmin / (vmax - vmin)

    print(f"Creating a diverging colormap centered at zero (position: {zero_norm:.2f})")
    
    # 定义色图的节点和颜色
    # 格式为: (位置比例, 颜色)
    return mcolors.LinearSegmentedColormap.from_list(
        "custom_blue_yellow_red",
        [(0, "blue"), (zero_norm, "yellow"), (1, "red")]
    )


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="Visualize a matrix as a 2D heatmap with a custom blue-yellow-red colormap.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-f', '--infile',
        required=True,
        metavar='matrix.dat',
        help="Path to the input matrix file (e.g., stiffness.dat)."
    )
    parser.add_argument(
        '-o', '--outfile',
        required=True,
        metavar='heatmap.png',
        help="Path to save the output heatmap image (e.g., heatmap.png)."
    )
    parser.add_argument(
        '--title',
        default='Matrix Heatmap',
        help="Title for the plot."
    )
    # 不再需要 --cmap 参数，因为我们总是使用自定义色图
    parser.add_argument(
        '--blocksize',
        type=int,
        default=None,
        help="Size of the blocks to draw major gridlines for (e.g., 6 for a 6x6 block)."
    )
    parser.add_argument(
        '--vmin',
        type=float,
        default=None,
        help="Minimum value for the color scale (optional, overrides auto-detection)."
    )
    parser.add_argument(
        '--vmax',
        type=float,
        default=None,
        help="Maximum value for the color scale (optional, overrides auto-detection)."
    )
    return parser.parse_args()

def main():
    """主执行函数"""
    args = parse_arguments()

    print(f"Reading matrix from: {args.infile}")
    try:
        matrix = np.loadtxt(args.infile)
    except FileNotFoundError:
        print(f"Error: Input file not found at '{args.infile}'")
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to read or parse the input file. Reason: {e}")
        sys.exit(1)

    print("Generating heatmap with custom colormap...")
    
    plt.figure(figsize=(10, 10))

    # --- 修改点: 使用自定义色图 ---
    # 1. 确定颜色范围
    vmin = args.vmin if args.vmin is not None else matrix.min()
    vmax = args.vmax if args.vmax is not None else matrix.max()
    
    # 2. 创建自定义色图对象
    custom_cmap = _create_custom_colormap(vmin, vmax)

    # 3. 将色图和范围传递给 seaborn
    ax = sns.heatmap(
        matrix,
        cmap=custom_cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.5,
        linecolor='white',
        cbar_kws={'orientation': 'horizontal', 'pad': 0.1, 'shrink': 0.7},
        square=True
    )
    
    ax.figure.axes[-1].set_xlabel('Value', fontsize=12)
    ax.set_title(args.title, fontsize=16, pad=20)
    
    dim = matrix.shape[0]

    if args.blocksize and args.blocksize > 0:
        num_blocks = (dim + args.blocksize - 1) // args.blocksize
        for i in range(1, num_blocks):
            pos = i * args.blocksize
            ax.axvline(pos, color='black', linewidth=2.0)
            ax.axhline(pos, color='black', linewidth=2.0)

        tick_locs = np.arange(dim) + 0.5
        tick_labels = [str((i % args.blocksize) + 1) for i in range(dim)]
        
        ax.set_xticks(tick_locs)
        ax.set_xticklabels(tick_labels, fontsize=9)
        ax.set_yticks(tick_locs)
        ax.set_yticklabels(tick_labels, fontsize=9, rotation=0)
        
        ax.tick_params(axis='both', which='both', length=0)
        ax.set_xlabel('')
        ax.set_ylabel('')
    else:
        ax.set_xlabel("Index", fontsize=12)
        ax.set_ylabel("Index", fontsize=12)

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    
    try:
        plt.savefig(args.outfile, dpi=300)
        print(f"Heatmap successfully saved to: {args.outfile}")
    except Exception as e:
        print(f"Error: Failed to save the figure. Reason: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
