import subprocess
import csv
import random
import re
import numpy as np

# 配置参数
NUM_TESTS = 1500
MAX_VECTOR_LEN = 1000
MAX_VALUE = 1000
OUTPUT_FILE = "results.csv"
TYPE = "int"

def generate_input_pair():
    """生成两行输入（每行一个向量）"""
    length = MAX_VECTOR_LEN
    if TYPE == "int":
        vec1 = [str(random.randint(-MAX_VALUE, MAX_VALUE)) for _ in range(length)]
        vec2 = [str(random.randint(-MAX_VALUE, MAX_VALUE)) for _ in range(length)]
        return f"{','.join(vec1)}\n{','.join(vec2)}"
    elif TYPE == "double":
        vec1 = [f"{random.uniform(-MAX_VALUE, MAX_VALUE):.10f}" for _ in range(length)]
        vec2 = [f"{random.uniform(-MAX_VALUE, MAX_VALUE):.10f}" for _ in range(length)]
        return f"{','.join(vec1)}\n{','.join(vec2)}"

def parse_program_output(output):
    """解析程序输出中的点积结果和耗时"""
    lines = output.split('\n')
    if len(lines) < 2:
        return None, None
    result_match = re.search(r'Dot product: (.*)', lines[0])
    time_match = re.search(r'Execution time: ([\d.]+) seconds', lines[1])
    if not result_match or not time_match:
        return None, None
    return result_match.group(1), float(time_match.group(1))

def run_test(input_data):
    """运行程序并返回输出"""
    proc = subprocess.run(
        ["./dotproductAos"],
        input=input_data,
        text=True,
        capture_output=True
    )
    return proc.stdout.strip()

def calculate_stats(times):
    """计算时间数据的平均值和方差"""
    if not times:
        return 0.0, 0.0
    times_array = np.array(times)
    mean = np.mean(times_array)
    variance = np.var(times_array, ddof=1)  # 无偏方差
    return mean, variance

def main():
    all_times = []  # 存储所有测试的耗时
    
    with open(OUTPUT_FILE, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["vector1", "vector2", "output", "time_seconds"])
        
        for i in range(NUM_TESTS):
            input_pair = generate_input_pair()
            lines = input_pair.split('\n')
            vec1, vec2 = lines[0], lines[1]
            
            try:
                output = run_test(input_pair)
                result, elapsed = parse_program_output(output)
                if result is None or elapsed is None:
                    print(f"解析失败: 输出=\n{output}")
                    continue
                
                # 记录数据
                writer.writerow([vec1, vec2, result, elapsed])
                all_times.append(elapsed)
                print(f"测试{i}完成: 耗时={elapsed:.9f} 秒")
                
            except Exception as e:
                print(f"测试失败: {e}")
    
    # 计算统计量
    mean, variance = calculate_stats(all_times)
    
    # 输出统计结果
    print("\n============= 统计结果 =============")
    print(f"测试次数: {len(all_times)}")
    print(f"平均耗时: {mean:.9f} 秒")
    print(f"耗时标准差: {np.sqrt(variance):.9f} 秒")


if __name__ == "__main__":
    main()