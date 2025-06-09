import cv2
import numpy as np
import time

# ================== 原始版本（使用cv2.add） ==================
def adjust_brightness_original(img, delta):
    """原始亮度调整方法"""
    start_time = time.perf_counter()
    result = cv2.add(img, delta)
    end_time = time.perf_counter()
    print(f"原始方法耗时：{(end_time - start_time)*1000:.2f}ms")
    return result

# ================== LUT优化版本 ==================
def adjust_brightness_lut(img, delta):
    """LUT优化亮度调整"""
    start_time = time.perf_counter()
    
    # 生成LUT表（核心优化点）
    lut = np.clip(np.arange(256, dtype=np.int16) + delta, 0, 255).astype(np.uint8)
    
    # 应用查找表
    result = cv2.LUT(img, lut)
    
    end_time = time.perf_counter()
    print(f"LUT方法耗时：{(end_time - start_time)*1000:.2f}ms")
    return result

# ================== 性能对比测试 ==================
def benchmark(image_path, delta, runs=100):
    """性能基准测试"""
    print(f"\n正在测试 {image_path} (运行次数：{runs})")
    
    # 加载测试图像
    img = cv2.imread(image_path)
    if img is None:
        print("无法加载图像")
        return

    # 预热运行（避免冷启动误差）
    _ = adjust_brightness_original(img.copy(), delta)
    _ = adjust_brightness_lut(img.copy(), delta)

    # 原始方法测试
    original_times = []
    for _ in range(runs):
        img_copy = img.copy()
        start = time.perf_counter()
        _ = adjust_brightness_original(img_copy, delta)
        original_times.append(time.perf_counter() - start)

    # LUT方法测试
    lut_times = []
    for _ in range(runs):
        img_copy = img.copy()
        start = time.perf_counter()
        _ = adjust_brightness_lut(img_copy, delta)
        lut_times.append(time.perf_counter() - start)

    # 统计结果
    avg_original = np.mean(original_times) * 1000
    avg_lut = np.mean(lut_times) * 1000
    improvement = (avg_original - avg_lut) / avg_original * 100

    print("\n测试结果：")
    print(f"原始方法平均耗时：{avg_original:.2f}ms")
    print(f"LUT方法平均耗时：{avg_lut:.2f}ms")
    print(f"性能提升：{improvement:.1f}%")

# ================== 使用示例 ==================
if __name__ == "__main__":
    # 测试参数设置
    test_image = "t4320.bmp"  # 替换为你的测试图像路径
    brightness_delta = 50          # 亮度调整值
    test_runs = 100                # 测试次数

    # 执行性能测试
    benchmark(test_image, brightness_delta, test_runs)

    # 可视化验证结果一致性
    img = cv2.imread(test_image)
    
    # 原始方法处理
    result_original = adjust_brightness_original(img.copy(), brightness_delta)
    
    # LUT方法处理
    result_lut = adjust_brightness_lut(img.copy(), brightness_delta)
    
    # 显示结果对比（按任意键切换）
    cv2.imshow("Original Result", result_original)
    cv2.waitKey(0)
    cv2.imshow("LUT Optimized Result", result_lut)
    cv2.waitKey(0)
    cv2.destroyAllWindows()