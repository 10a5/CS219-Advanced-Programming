import subprocess
import csv
import random
import re
import numpy as np

# 配置参数（测试时先用小数据量）
NUM_TESTS = 1500
MAX_VECTOR_LEN = 1000 
MAX_VALUE = 1000
OUTPUT_CSV = "results.csv"
DATA_TYPE = "double"

def generate_input_pair():
    """生成两行输入数据（用空格分隔）"""
    length = MAX_VECTOR_LEN
    if DATA_TYPE == "int":
        vec1 = [str(random.randint(-MAX_VALUE, MAX_VALUE)) for _ in range(length)]
        vec2 = [str(random.randint(-MAX_VALUE, MAX_VALUE)) for _ in range(length)]
    else:
        vec1 = [f"{random.uniform(-MAX_VALUE, MAX_VALUE):.10f}" for _ in range(length)]
        vec2 = [f"{random.uniform(-MAX_VALUE, MAX_VALUE):.10f}" for _ in range(length)]
    return f"{' '.join(vec1)}\n{' '.join(vec2)}"

def parse_program_output(output):
    """解析程序输出（兼容不同换行符）"""
    result_match = re.search(r'Dot product:\s*(.*?)\s*\n', output)
    time_match = re.search(r'Time elapsed:\s*([\d.]+)\s*seconds', output)
    return (
        result_match.group(1).strip() if result_match else None,
        float(time_match.group(1)) if time_match else None
    )

def run_test(input_data):
    """运行Java程序（分配足够内存）"""
    proc = subprocess.run(
        ["java", "-Xms512m", "-Xmx2048m", "Dotproduct"],
        input=input_data,
        text=True,
        capture_output=True
    )
    if proc.returncode != 0:
        print(f"Java程序错误: {proc.stderr.strip()}")
        return ""
    return proc.stdout.strip()

def calculate_stats(times):
    """计算统计量"""
    times_array = np.array(times)
    return np.mean(times_array), np.var(times_array, ddof=1)

def main():
    all_times = []
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Vector1", "Vector2", "Result", "Time(seconds)"])
        
        for _ in range(NUM_TESTS):
            input_data = generate_input_pair()
            vec1, vec2 = input_data.split("\n")
            
            try:
                output = run_test(input_data)
                result, elapsed = parse_program_output(output)
                if not result or not elapsed:
                    print(f"解析失败: {output}")
                    continue
                
                writer.writerow([vec1, vec2, result, elapsed])
                all_times.append(elapsed)
                # print(f"测试完成: 耗时 {elapsed:.9f} 秒")
            except Exception as e:
                print(f"测试失败: {str(e)}")
    
    if all_times:
        mean, var = calculate_stats(all_times)
        print("\n======== 统计结果 ========")
        print(f"总测试次数: {len(all_times)}")
        print(f"平均耗时  : {mean:.9f} 秒")
        print(f"标准差    : {np.sqrt(var):.9f} 秒")

if __name__ == "__main__":
    main()