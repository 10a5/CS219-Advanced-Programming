import matplotlib.pyplot as plt
import numpy as np
import os
os.chdir(os.path.dirname(__file__))

# === 数据配置 ===
sizes = [
    '1000×600', 
    '2000×1200', 
    '5000×3000',
    '10000×6000',
    '15000×9000',
    '20000×12000'
]

before = [0.07, 2.05, 11.45, 42.46, 92.14, 158.94]
after = [0.06, 0.22, 4.07, 15.65, 32.2, 56.72]

x_values = [int(w)/1e3 for size in sizes for w,h in [size.split('×')]]

# === 可视化配置 ===
plt.figure(figsize=(12,7))

lines = [
    ('before', before, 'blue', '-', 'o'),
    ('after', after, 'purple', '--', 's')
]

for label, data, color, linestyle, marker in lines:
    plt.plot(x_values, data, 
             label=label, 
             color=color,
             linestyle=linestyle,
             marker=marker,
             markersize=8,
             linewidth=2,
             alpha=0.8)

# === 坐标轴优化 ===
plt.title('Processing Time with Uniform Y-Axis', fontsize=14, pad=20)
plt.xlabel('Image Size (Megapixels)', fontsize=12)
plt.ylabel('Processing Time (ms)', fontsize=12)
plt.xticks(x_values, sizes, rotation=35, ha='right')

# y轴均匀分布设置
ax = plt.gca()
ax.yaxis.set_major_locator(plt.MultipleLocator(100))  # 主刻度每100ms
ax.yaxis.set_minor_locator(plt.MultipleLocator(50))   # 次刻度每50ms
plt.ylim(0, 200)  # 固定显示范围

plt.grid(True, linestyle='--', alpha=0.6, which='both')  # 显示主次网格
plt.legend(frameon=True, fontsize=10)

# === 输出 ===
plt.tight_layout()
plt.savefig('uniform_y_axis2.png', dpi=300, bbox_inches='tight')