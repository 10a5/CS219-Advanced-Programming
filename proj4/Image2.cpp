// Image.cpp
#include "One.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace SimpleImage {

// --- 私有方法实现 ---
void Image::releaseResource() {
    if (m_refCount) {
        (*m_refCount)--;
        if (*m_refCount == 0) {
            delete[] m_data;
            delete m_refCount;
        }
        m_data = nullptr;
        m_refCount = nullptr;
        m_width = m_height = m_channels = 0;
    }
}

// void Image::detach() {
//     if (m_refCount && *m_refCount > 1) {
//         unsigned char* newData = new unsigned char[m_width * m_height * m_channels];
//         std::memcpy(newData, m_data, m_width * m_height * m_channels);
//         releaseResource();
//         m_data = newData;
//         m_refCount = new int(1);
//     }
// }

void Image::checkBounds(int x, int y, int channel) const {
    if (x < 0 || x >= m_width || y < 0 || y >= m_height || channel < 0 || channel >= m_channels) {
        throw std::out_of_range("Pixel access out of bounds");
    }
}

// --- 构造函数与析构函数 ---
Image::Image() = default;

Image::Image(int width, int height, int channels) 
    : m_width(width), m_height(height), m_channels(channels) {
    if (width <= 0 || height <= 0 || channels != 3 && channels != 1) {
        throw std::invalid_argument("Invalid image dimensions when creating Image");
    }
    m_data = new unsigned char[width * height * channels];
    m_refCount = new int(1);
}

Image::Image(const Image& other) 
    : m_width(other.m_width), m_height(other.m_height), m_channels(other.m_channels),
      m_data(other.m_data), m_refCount(other.m_refCount) {
    if (m_refCount) (*m_refCount)++;
}

Image::Image(Image&& other) 
    : m_width(other.m_width), m_height(other.m_height), m_channels(other.m_channels),
      m_data(other.m_data), m_refCount(other.m_refCount) {
    other.m_width = other.m_height = other.m_channels = 0;
    other.m_data = nullptr;
    other.m_refCount = nullptr;
}

Image::~Image() {
    releaseResource();
}

// --- 赋值运算符 ---
Image& Image::operator=(const Image& other) {
    if (this != &other) {
        releaseResource();
        m_width = other.m_width;
        m_height = other.m_height;
        m_channels = other.m_channels;
        m_data = other.m_data;
        m_refCount = other.m_refCount;
        if (m_refCount) (*m_refCount)++;
    }
    return *this;
}

Image& Image::operator=(Image&& other) {
    if (this != &other) {
        releaseResource();
        m_width = other.m_width;
        m_height = other.m_height;
        m_channels = other.m_channels;
        m_data = other.m_data;
        m_refCount = other.m_refCount;
        other.m_width = other.m_height = other.m_channels = 0;
        other.m_data = nullptr;
        other.m_refCount = nullptr;
    }
    return *this;
}

// --- 图像加载/保存 ---
bool Image::load(const std::string& filename) {
    releaseResource();
    
    int channels = 0;
    unsigned char* stb_data = stbi_load(filename.c_str(), &m_width, &m_height, &channels, 0);
    if (!stb_data) {
        throw std::invalid_argument ("Failed to load image: " + filename);
        return false;
    }
    
    m_channels = channels;
    size_t data_size = m_width * m_height * m_channels;
    m_data = new unsigned char[data_size];
    std::memcpy(m_data, stb_data, data_size);
    stbi_image_free(stb_data);
    
    m_refCount = new int(1);
    return true;
}

bool Image::save(const std::string& filename) const {
    if (!m_data) return false;
    
    std::string ext = filename.substr(filename.find_last_of(".") + 1);
    int success = 0;
    
    if (ext == "png") {
        success = stbi_write_png(filename.c_str(), m_width, m_height, m_channels, m_data, m_width * m_channels);
    } else if (ext == "jpg" || ext == "jpeg") {
        success = stbi_write_jpg(filename.c_str(), m_width, m_height, m_channels, m_data, 90);
    } else if (ext == "bmp") {
        success = stbi_write_bmp(filename.c_str(), m_width, m_height, m_channels, m_data);
    } else {
        throw std::invalid_argument("Unsupported image format: " + ext);
        return false;
    }
    
    if (!success) {
        throw std::invalid_argument("Failed to save image: " + filename);
        return false;
    }
    return true;
}

// --- 图像处理操作 ---
void Image::adjustBrightness(int value) {
    if (!m_data) return;
    // detach();
    
    #pragma omp parallel for
    for (int i = 0; i < m_width * m_height * m_channels; ++i) {
        int new_val = static_cast<int>(m_data[i]) + value;
        m_data[i] = static_cast<unsigned char>(std::clamp(new_val, 0, 255));
    }
}

void Image::blend(const Image& other, float alpha) {
    if (!m_data || !other.m_data) return;
    if (m_width != other.m_width || m_height != other.m_height || m_channels != other.m_channels) {
        throw std::invalid_argument("Images must have the same dimensions");
    }
    
    // detach();
    float beta = 1.0f - alpha;
    
    #pragma omp parallel for
    for (int i = 0; i < m_width * m_height * m_channels; ++i) {
        float blended = alpha * m_data[i] + beta * other.m_data[i];
        m_data[i] = static_cast<unsigned char>(std::clamp(blended, 0.0f, 255.0f));
    }
}

// --- 像素访问 ---
unsigned char& Image::at(int x, int y, int channel) {
    checkBounds(x, y, channel);
    // detach();
    return m_data[(y * m_width + x) * m_channels + channel];
}

const unsigned char& Image::at(int x, int y, int channel) const {
    checkBounds(x, y, channel);
    return m_data[(y * m_width + x) * m_channels + channel];
}

} // namespace SimpleImage