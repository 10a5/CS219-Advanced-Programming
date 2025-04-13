import matplotlib.pyplot as plt
import numpy as np
import os
os.chdir(os.path.dirname(__file__))

# === 数据配置 ===
sizes = [
    '1000×600', 
    '2000×1200', 
    '3000×1800',
    '4000×2400',
    '5000×3000',
    '8000×4800',
    '10000×6000',
    '12000×7200',
    '15000×9000',
    '20000×12000'
]

interlaced = [0.93, 3.95, 9.36, 17.32, 23.92, 60.50, 94.36, 134.26, 210.24, 374.45]
planar = [1.49, 6.42, 13.65, 24.19, 38.03, 96.14, 150.28, 217.63, 344.31, 607.79]
opt_interlaced = [0.75, 2.94, 7.55, 13.70, 21.61, 56.48, 86.46, 125.58, 197.37, 351.69]
opt_planar = [0.78, 2.97, 6.88, 13.03, 20.98, 55.52, 84.99, 123.39, 196.27, 349.25]

x_values = [int(w)/1e3 for size in sizes for w,h in [size.split('×')]]

# === 可视化配置 ===
plt.figure(figsize=(12,7))

lines = [
    ('Interlaced', interlaced, 'blue', '-', 'o'),
    ('Planar', planar, 'green', '-', 's'),
    ('Optimized Interlaced', opt_interlaced, 'red', '--', 'o'),
    ('Optimized Planar', opt_planar, 'purple', '--', 's')
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
plt.ylim(0, 650)  # 固定显示范围

plt.grid(True, linestyle='--', alpha=0.6, which='both')  # 显示主次网格
plt.legend(frameon=True, fontsize=10)

# === 输出 ===
plt.tight_layout()
plt.savefig('uniform_y_axis.png', dpi=300, bbox_inches='tight')