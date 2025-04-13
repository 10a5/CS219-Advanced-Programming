import os
import argparse
from PIL import Image

def convert_image(input_path, output_path=None):
    """
    支持多种格式互转：JPG/PNG -> BMP，JPG -> PNG
    :param input_path: 输入文件路径（或目录）
    :param output_path: 输出文件路径（或目录），默认为输入路径的同名文件
    """
    # 检查输入路径类型（文件或目录）
    if os.path.isdir(input_path):
        input_files = [
            f for f in os.listdir(input_path)
            if f.lower().endswith((".png", ".jpg", ".jpeg"))
        ]
        if not input_files:
            raise ValueError("目录中没有支持的图像文件（仅支持 PNG/JPG）")

        # 自动创建输出目录
        output_dir = output_path if output_path else input_path
        os.makedirs(output_dir, exist_ok=True)

        for file in input_files:
            input_file = os.path.join(input_path, file)
            base_name = os.path.splitext(file)[0]
            
            # 自动推断输出格式
            if output_path and os.path.isdir(output_path):
                output_file = os.path.join(output_dir, base_name + _get_output_ext(input_file, None))
            else:
                output_file = os.path.join(output_dir, base_name + _get_output_ext(input_file, output_path))
            
            _convert_single_file(input_file, output_file)
        
        print(f"转换完成！共转换 {len(input_files)} 个文件")
    else:
        # 单文件转换
        if not output_path:
            output_path = os.path.splitext(input_path)[0] + _get_output_ext(input_path, None)
        
        _convert_single_file(input_path, output_path)
        print(f"转换完成！输出文件: {output_path}")

def _get_output_ext(input_path, output_path):
    """根据输入和输出路径推断目标格式"""
    if output_path:
        ext = os.path.splitext(output_path)[1].lower()
        if ext not in (".bmp", ".png", ".jpg", ".jpeg"):
            raise ValueError("输出格式不支持（仅支持 .bmp/.png/.jpg）")
        return ext
    else:
        # 默认转 BMP
        return ".bmp"

def _convert_single_file(input_file, output_file):
    """处理单个文件转换"""
    try:
        with Image.open(input_file) as img:
            # 处理透明度（BMP/JPG 不支持透明通道）
            output_ext = os.path.splitext(output_file)[1].lower()
            if img.mode in ("RGBA", "LA") and output_ext in (".bmp", ".jpg", ".jpeg"):
                background = Image.new("RGB", img.size, (255, 255, 255))  # 透明背景转为白色
                background.paste(img, mask=img.split()[-1])
                img = background
            
            # 保存为对应格式
            if output_ext == ".bmp":
                img.save(output_file, "BMP")
            elif output_ext in (".jpg", ".jpeg"):
                img.save(output_file, "JPEG", quality=95)  # 保持高质量
            elif output_ext == ".png":
                img.save(output_file, "PNG")
            else:
                raise ValueError("不支持的输出格式")
    except Exception as e:
        raise RuntimeError(f"文件转换失败: {str(e)}")

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    parser = argparse.ArgumentParser(description="图像格式转换工具（支持 PNG/JPG/BMP）")
    parser.add_argument("-i", "--input", required=True, help="输入文件或目录路径")
    parser.add_argument("-o", "--output", help="输出文件或目录路径")
    args = parser.parse_args()

    try:
        convert_image(args.input, args.output)
    except Exception as e:
        print(f"错误: {str(e)}")
        print("使用示例:")
        print("  JPG 转 BMP: python convert.py -i input.jpg -o output.bmp")
        print("  PNG 转 JPG: python convert.py -i input.png -o output.jpg")
        print("  JPG 转 PNG: python convert.py -i input.jpg -o output.png")
        print("  批量转换:    python convert.py -i ./input_images/ -o ./output_folder/")