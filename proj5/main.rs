use clap::Parser;
use image::{ImageBuffer, RgbImage};
use rayon::prelude::*;
use std::time::Instant;

#[derive(Parser, Debug, Clone)]
#[clap(version = "1.0", author = "Your Name")]
struct Args {
    /// 输入图像路径
    #[clap(short = 'i', long = "input", required = true)]
    input_file: String,

    /// 输出图像路径
    #[clap(short = 'o', long = "output", required = true)]
    output_file: String,

    /// 亮度调节值（-255到255）
    #[clap(short = 'b', long = "brightness", required = true)]
    brightness_value: i32,

    /// 显示运行时间
    #[clap(short = 't', long = "time")]
    show_time: bool,
    
    /// 使用串行处理（默认并行）
    #[clap(long = "serial")]
    use_serial: bool,
    
    /// 运行性能测试（比较串行和并行）
    #[clap(long = "benchmark")]
    run_benchmark: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let start_time = Instant::now();

    // 读取输入图像
    let img = image::open(&args.input_file)?.to_rgb8();
    let (width, height) = img.dimensions();
    println!("图像尺寸: {}x{} ({} 像素)", width, height, width * height);
    
    let result = handle_brightness(&args, img);
    
    if args.show_time {
        let total_duration = start_time.elapsed();
        println!("\n总运行时间: {:.6}秒", total_duration.as_secs_f64());
    }

    result
}

fn handle_brightness(args: &Args, img: RgbImage) -> Result<(), Box<dyn std::error::Error>> {
    let brightness = args.brightness_value;
    
    if args.run_benchmark {
        // 性能对比测试
        benchmark_brightness(&img, brightness);
    } else {
        // 正常处理
        let start_time = Instant::now();
        
        let adjusted = if args.use_serial {
            adjust_brightness_serial(&img, brightness)
        } else {
            adjust_brightness_parallel(&img, brightness)
        };
        
        let duration = start_time.elapsed();
        
        adjusted.save(&args.output_file)?;
        println!("结果已保存至: {}", args.output_file);
        
        if args.show_time {
            let mode = if args.use_serial { "串行" } else { "并行" };
            println!("亮度调节处理时间({}模式): {:.6}秒", mode, duration.as_secs_f64());
        }
    }
    
    Ok(())
}

/// 性能对比测试
fn benchmark_brightness(img: &RgbImage, brightness: i32) {
    println!("\n开始亮度调节性能测试...");
    
    // 预热（确保CPU达到最大频率）
    let _warmup = adjust_brightness_serial(img, brightness);
    
    // 串行处理测试
    let serial_start = Instant::now();
    let _serial_result = adjust_brightness_serial(img, brightness);
    let serial_duration = serial_start.elapsed();
    
    // 并行处理测试
    let parallel_start = Instant::now();
    let _parallel_result = adjust_brightness_parallel(img, brightness);
    let parallel_duration = parallel_start.elapsed();
    
    // 输出结果
    println!("串行处理时间: {:.6}秒", serial_duration.as_secs_f64());
    println!("并行处理时间: {:.6}秒", parallel_duration.as_secs_f64());
    println!("加速比: {:.2}x", serial_duration.as_secs_f64() / parallel_duration.as_secs_f64());
    
    // 验证结果一致性
    let serial_pixels = _serial_result.pixels().count();
    let parallel_pixels = _parallel_result.pixels().count();
    println!("处理像素数验证: 串行={}, 并行={}", serial_pixels, parallel_pixels);
}

/// 串行亮度调节
fn adjust_brightness_serial(img: &RgbImage, value: i32) -> RgbImage {
    let (width, height) = img.dimensions();
    let delta = value.clamp(-128, 127) as i8;
    
    let mut output = ImageBuffer::new(width, height);
    
    for y in 0..height {
        for x in 0..width {
            let pixel = img.get_pixel(x, y);
            let mut channels = pixel.0;
            
            channels[0] = channels[0].saturating_add_signed(delta);
            channels[1] = channels[1].saturating_add_signed(delta);
            channels[2] = channels[2].saturating_add_signed(delta);
            
            *output.get_pixel_mut(x, y) = image::Rgb(channels);
        }
    }
    
    output
}

/// 并行亮度调节
fn adjust_brightness_parallel(img: &RgbImage, value: i32) -> RgbImage {
    let (width, height) = img.dimensions();
    let delta = value.clamp(-128, 127) as i8;
    
    // 创建输出缓冲区
    let mut output = ImageBuffer::new(width, height);
    
    // 使用并行迭代器处理像素
    output.par_chunks_exact_mut(3).enumerate().for_each(|(i, pixel)| {
        let x = (i % width as usize) as u32;
        let y = (i / width as usize) as u32;
        
        let src_pixel = img.get_pixel(x, y);
        pixel[0] = src_pixel[0].saturating_add_signed(delta);
        pixel[1] = src_pixel[1].saturating_add_signed(delta);
        pixel[2] = src_pixel[2].saturating_add_signed(delta);
    });
    
    output
}