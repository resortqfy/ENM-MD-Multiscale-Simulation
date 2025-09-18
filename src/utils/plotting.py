import matplotlib.pyplot as plt
import seaborn as sns
import py3Dmol  # 需额外安装

def plot_bfactors(bfactors: np.ndarray, output_file: str = None):
    """绘制B因子"""
    sns.lineplot(x=range(len(bfactors)), y=bfactors)
    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()

def visualize_modes(coords: np.ndarray, modes: np.ndarray, index: int = 0):
    """可视化ENM模式"""
    # 简化，使用py3Dmol显示
    view = py3Dmol.view()
    # 添加模型和动画逻辑
    view.show()
