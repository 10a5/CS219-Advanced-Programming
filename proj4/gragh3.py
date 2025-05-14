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

# 使用表格中提供的新时间数据
before = [0.92, 4.05, 21.51, 83.07, 143.55, 329.75]  # 优化前时间
after = [0.21, 1.01, 4.66, 17.27, 38.01, 66.40]      # 优化后时间

# 计算x轴值（单位：百万像素）
x_values = [int(w)/1e3 for size in sizes for w,h in [size.split('×')]]

# === 可视化配置 ===
plt.figure(figsize=(12,7))

# 定义线条样式
lines = [
    ('Before Optimization', before, '#1f77b4', '-', 'o'),  # 蓝色
    ('After Optimization', after, '#ff7f0e', '--', 's')    # 橙色
]

# 绘制折线图
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
plt.title('Image Processing Time Comparison', fontsize=14, pad=20)
plt.xlabel('Image Size (Megapixels)', fontsize=12)
plt.ylabel('Processing Time (ms)', fontsize=12)
plt.xticks(x_values, sizes, rotation=35, ha='right')

# 调整y轴范围
ax = plt.gca()
ax.set_ylim(0, 350)  # 根据最大时间329.75调整
ax.yaxis.set_major_locator(plt.MultipleLocator(50))  # 主刻度每50ms
ax.yaxis.set_minor_locator(plt.MultipleLocator(25))  # 次刻度每25ms

# 添加网格
plt.grid(True, linestyle='--', alpha=0.6, which='both')

# 显示图例
plt.legend(frameon=True, fontsize=10, loc='upper left')

# === 输出 ===
plt.tight_layout()
plt.savefig('processing_time_comparison1.png', dpi=300, bbox_inches='tight')
plt.show()