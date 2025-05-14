#include "One.h"
#include <iostream>
#include <cassert>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;
using namespace SimpleImage;

// 辅助函数：创建测试文件
void create_test_file(const std::string& filename, int size = 4) {
    std::ofstream file(filename, std::ios::binary);
    unsigned char data[] = {255,0,0,0,255,0,0,0,255}; // 3x3 RGB
    file.write(reinterpret_cast<char*>(data), sizeof(data));
}

// 测试基类Image功能
void test_base_image() {
    // 默认构造
    Image img1;
    assert(img1.width() == 0 && img1.height() == 0);

    // 参数构造
    Image img2(100, 200, 3);
    assert(img2.width() == 100 && img2.channels() == 3);

    // 拷贝构造
    Image img3 = img2;
    assert(img3.data() != img2.data()); // 深拷贝检查

    // 赋值操作符
    Image img4;
    img4 = img3;
    assert(img4.height() == 200 && img4.data() != img3.data());
}

// 测试JpegImage
void test_jpeg_image() {
    create_test_file("test.jpg");
    JpegImage jpg1("test.jpg");
    assert(jpg1.width() == 3 && jpg1.channels() == 3);
    assert(jpg1.save("output.jpg"));

    JpegImage jpg2(800, 600);
    jpg2.setQuality(50);
    assert(jpg2.save("empty.jpg"));

    bool exception_thrown = false;
    try {
        JpegImage jpg3("nonexist.jpg");
    } catch (...) {
        exception_thrown = true;
    }
    assert(exception_thrown);
}

// 测试GrayscaleImage
void test_grayscale_image() {
    create_test_file("test_gray.jpg");
    GrayscaleImage gray1("test_gray.jpg");
    assert(gray1.channels() == 1 && gray1.getPixel(0,0) == 76); // 0.299*255 ≈76

    GrayscaleImage gray2(10, 10);
    gray2.setPixel(5,5, 128);
    assert(gray2.getPixel(5,5) == 128);
}

// 测试PngImage
void test_png_image() {
    bool exception_thrown = false;
    try {
        PngImage png1(100, 100, 2); // 无效通道
    } catch (...) {
        exception_thrown = true;
    }
    assert(exception_thrown);

    PngImage png2(256, 256, 4);
    assert(png2.save("test.png")); // 需要Image的save实现
}

// 测试Image256
void test_image256() {
    Image256 img1(100, 100, 1);
    img1.adjustBrightness(50);
    assert(img1.data8()[0] <= 255);

    create_test_file("test_256.png");
    Image256 img2("test_256.png");
    assert(img2.channels() == 3); // STB加载的通道数可能为3

    Image256 img3(50,50,3);
    img3.blend(img2, 0.7f);
}

int main() {
    test_base_image();
    test_jpeg_image();
    test_grayscale_image();
    test_png_image();
    test_image256();

    // 清理测试文件
    fs::remove("test.jpg");
    fs::remove("output.jpg");
    fs::remove("empty.jpg");
    fs::remove("test_gray.jpg");
    fs::remove("test.png");
    fs::remove("test_256.png");

    std::cout << "All tests passed!" << std::endl;
    return 0;
}