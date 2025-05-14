#include <vector>
#include <cstdint>
#include <string>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// 生成单通道（灰度）测试图像
void create_grayscale_test_image(const std::string& filename, 
                                int width, 
                                int height, 
                                uint8_t pattern_type = 0) 
{
    std::vector<uint8_t> pixels(width * height);

    // 生成不同测试模式
    switch(pattern_type) {
        case 0: // 从左到右线性渐变（黑到白）
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    pixels[y * width + x] = static_cast<uint8_t>(255 * x / (width-1));
                }
            }
            break;

        case 1: // 棋盘格模式
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    bool checker = ((x/32 + y/32) % 2) == 0;
                    pixels[y * width + x] = checker ? 255 : 0;
                }
            }
            break;

        case 2: // 中心渐变（白到黑）
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    float dx = (x - width/2.0f) / (width/2.0f);
                    float dy = (y - height/2.0f) / (height/2.0f);
                    float dist = sqrt(dx*dx + dy*dy);
                    pixels[y * width + x] = static_cast<uint8_t>(255 * (1.0f - dist));
                }
            }
            break;

        default: // 全白背景
            std::fill(pixels.begin(), pixels.end(), 255);
    }

    // 保存为PNG（单通道）
    stbi_write_png(filename.c_str(), 
                  width, 
                  height, 
                  1,  // 单通道
                  pixels.data(), 
                  width); // 每行字节数=width*1
}

int main() {
    // 生成三种测试图
    create_grayscale_test_image("gradient.png", 512, 512, 0);
    create_grayscale_test_image("checkerboard.png", 512, 512, 1);
    create_grayscale_test_image("radial_gradient.png", 512, 512, 2);
    return 0;
}