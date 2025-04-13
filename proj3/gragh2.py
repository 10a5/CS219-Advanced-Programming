import matplotlib.pyplot as plt
import numpy as np
import os
os.chdir(os.path.dirname(__file__))
# 数据准备
operations = [
    'adjust brightness', 'average', 'resize', 'contrast', 'grayscale',
    'invert', 'crop', 'rotate_90_degrees', 'rotate_180_degrees',
    'vignette', 'rgb_shift'
]
pre_optimization = [
    1368.77, 196.59, 186.55, 535.65, 71.09,
    12.22, 213.13, 483.60, 96.18, 350.25, 171.69
]
post_optimization = [
    9.53, 38.58, 47.34, 92.6, 8.81,
    7.59, 30.15, 346.12, 23.75, 8.06, 12.34
]

# 创建位置数组
x = np.arange(len(operations))
width = 0.35  # 柱状图宽度

# 创建图形
plt.figure(figsize=(14, 8))

# 绘制柱状图
bars1 = plt.bar(x - width/2, pre_optimization, width, label='before optimization')
bars2 = plt.bar(x + width/2, post_optimization, width, label='after optimization')

# 添加标签、标题和图例
plt.ylabel('running time(ms)')
plt.title('Running Time of Image Processing Operations')
plt.xlabel('Image Processing Operations')
plt.xticks(x, operations, rotation=45, ha='right')
plt.legend()

# 设置对数坐标轴（因时间差异较大）
plt.yscale('log')
plt.ylim(1, 10000)
plt.yticks([1, 10, 100, 1000, 10000], ['1', '10', '100', '1000', '10000'])

# 添加数值标签
def add_labels(bars):
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{height:.1f}',
                 ha='center', va='bottom',
                 fontsize=8)

add_labels(bars1)
add_labels(bars2)

# 调整布局
plt.tight_layout()
plt.grid(axis='y', alpha=0.5, linestyle='--')

# 显示图形
plt.show()
plt.savefig('gragh2', dpi=300, bbox_inches='tight')
