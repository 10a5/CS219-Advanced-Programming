#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

#pragma pack(push, 1) // 确保结构体按 1 字节对齐

struct BMPFileHeader {
    uint16_t bfType = 0x4D42; // "BM"
    uint32_t bfSize;
    uint16_t bfReserved1 = 0;
    uint16_t bfReserved2 = 0;
    uint32_t bfOffBits = 54; // 头部大小 14 + 40
};

struct BMPInfoHeader {
    uint32_t biSize = 40;
    int32_t biWidth;
    int32_t biHeight;
    uint16_t biPlanes = 1;
    uint16_t biBitCount = 24; // 24 位颜色深度
    uint32_t biCompression = 0;
    uint32_t biSizeImage;
    int32_t biXPelsPerMeter = 0;
    int32_t biYPelsPerMeter = 0;
    uint32_t biClrUsed = 0;
    uint32_t biClrImportant = 0;
};

#pragma pack(pop) // 恢复默认对齐

void generateBMP(const char* filename, int width, int height) {
    BMPFileHeader fileHeader;
    BMPInfoHeader infoHeader;
    
    infoHeader.biWidth = width;
    infoHeader.biHeight = height;
    
    int rowSize = (width * 3 + 3) & ~3; // 4 字节对齐
    infoHeader.biSizeImage = rowSize * height;
    fileHeader.bfSize = fileHeader.bfOffBits + infoHeader.biSizeImage;
    
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "无法创建 BMP 文件" << std::endl;
        return;
    }
    
    file.write(reinterpret_cast<char*>(&fileHeader), sizeof(fileHeader));
    file.write(reinterpret_cast<char*>(&infoHeader), sizeof(infoHeader));
    
    std::vector<uint8_t> row(rowSize, 0);
    for (int i = 0; i < width * 3; i += 3) {
        row[i] = 0;     // 蓝色分量
        row[i + 1] = 0; // 绿色分量
        row[i + 2] = 255; // 红色分量
    }
    
    for (int y = 0; y < height; y++) {
        file.write(reinterpret_cast<char*>(row.data()), rowSize);
    }
    
    file.close();
    std::cout << "BMP 文件生成成功: " << filename << std::endl;
}

int main() {
    generateBMP("output.bmp", 100, 100);
    return 0;
}
