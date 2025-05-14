// Image.cpp
#include "One.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <immintrin.h>

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

// Image::Image(Image&& other) 
//     : m_width(other.m_width), m_height(other.m_height), m_channels(other.m_channels),
//       m_data(other.m_data), m_refCount(other.m_refCount) {
//     other.m_width = other.m_height = other.m_channels = 0;
//     other.m_data = nullptr;
//     other.m_refCount = nullptr;
// }

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

// Image& Image::operator=(Image&& other) {
//     if (this != &other) {
//         releaseResource();
//         m_width = other.m_width;
//         m_height = other.m_height;
//         m_channels = other.m_channels;
//         m_data = other.m_data;
//         m_refCount = other.m_refCount;
//         other.m_width = other.m_height = other.m_channels = 0;
//         other.m_data = nullptr;
//         other.m_refCount = nullptr;
//     }
//     return *this;
// }

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
    if(this->m_refCount == nullptr) {
        throw std::invalid_argument("Uninitiallized image!");
    }
    if (!m_data || value == 0) return;
    // detach();  // 确保数据独立（取消注释）

    const int total = m_width * m_height * m_channels;
    constexpr int simd_step = 32;  // 每次处理 32 字节（256 位）

    // 将亮度值限制在 [-255, 255] 并转为 16 位有符号
    value = std::clamp(value, -255, 255);
    const __m256i v_value = _mm256_set1_epi16(static_cast<short>(value));

    #pragma omp parallel for
    for (int i = 0; i <= total - simd_step; i += simd_step) {
        // 加载 32 字节数据
        __m256i data = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(m_data + i));

        // 将 8 位无符号像素扩展为 16 位有符号（分高低两部分）
        __m256i lo = _mm256_cvtepu8_epi16(_mm256_castsi256_si128(data));  // 低 128 位
        __m256i hi = _mm256_cvtepu8_epi16(_mm256_extracti128_si256(data, 1));  // 高 128 位

        // 执行 16 位无符号饱和加法（自动限制在 0~255）
        lo = _mm256_adds_epu16(lo, v_value);
        hi = _mm256_adds_epu16(hi, v_value);

        // 将 16 位结果压缩回 8 位无符号
        __m256i blended = _mm256_packus_epi16(lo, hi);
        blended = _mm256_permute4x64_epi64(blended, 0xD8);  // 修复通道顺序

        // 存储结果
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(m_data + i), blended);
    }

    // 处理剩余像素（标量循环）
    int remainder = total % simd_step;
    int tail_start = total - remainder;
    for (int i = tail_start; i < total; ++i) {
        int new_val = static_cast<int>(m_data[i]) + value;
        m_data[i] = static_cast<unsigned char>(std::clamp(new_val, 0, 255));
    }
}

void Image::blend(const Image& other, float alpha) {
    if(this->m_refCount == nullptr || other.m_refCount == nullptr) {
        throw std::invalid_argument("Uninitiallized image!");
    }
    if (!m_data || !other.m_data) return;
    if (m_width != other.m_width || m_height != other.m_height || m_channels != other.m_channels) {
        throw std::invalid_argument("Images must have the same dimensions");
    }
    
    // detach();  // 若启用写时复制，取消注释
    const float beta = 1.0f - alpha;
    const int total = m_width * m_height * m_channels;
    constexpr int simd_step = 32;  // 每次处理 32 字节（32像素 for 8UC1, 8像素 for 8UC4）

    // 预计算 SIMD 常量
    const __m256 alpha_vec = _mm256_set1_ps(alpha);
    const __m256 beta_vec = _mm256_set1_ps(beta);
    const __m256 zero = _mm256_setzero_ps();
    const __m256 max_val = _mm256_set1_ps(255.0f);

    // 主循环：每次处理 32 字节（SIMD 友好）
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < total - simd_step + 1; i += simd_step) {
        // 加载当前图像和另一图像的 32 字节数据（未对齐）
        __m256i this_data = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(m_data + i));
        __m256i other_data = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(other.m_data + i));

        // 分离高低 16 字节并处理
        for (int offset = 0; offset < 32; offset += 16) {
            // 提取当前图像的 16 字节（128-bit）
            __m128i this_chunk = _mm256_extracti128_si256(this_data, offset / 16);
            // 转换为 8 个 32-bit 整数并转浮点
            __m256i this_lo = _mm256_cvtepu8_epi32(this_chunk);
            __m256 this_lo_f = _mm256_cvtepi32_ps(this_lo);
            __m256i this_hi = _mm256_cvtepu8_epi32(_mm_unpackhi_epi64(this_chunk, this_chunk));
            __m256 this_hi_f = _mm256_cvtepi32_ps(this_hi);

            // 提取另一图像的 16 字节并同样处理
            __m128i other_chunk = _mm256_extracti128_si256(other_data, offset / 16);
            __m256i other_lo = _mm256_cvtepu8_epi32(other_chunk);
            __m256 other_lo_f = _mm256_cvtepi32_ps(other_lo);
            __m256i other_hi = _mm256_cvtepu8_epi32(_mm_unpackhi_epi64(other_chunk, other_chunk));
            __m256 other_hi_f = _mm256_cvtepi32_ps(other_hi);

            // 混合计算：alpha * this + beta * other
            __m256 blended_lo = _mm256_fmadd_ps(alpha_vec, this_lo_f, _mm256_mul_ps(beta_vec, other_lo_f));
            __m256 blended_hi = _mm256_fmadd_ps(alpha_vec, this_hi_f, _mm256_mul_ps(beta_vec, other_hi_f));

            // 截断到 [0, 255]
            blended_lo = _mm256_min_ps(_mm256_max_ps(blended_lo, zero), max_val);
            blended_hi = _mm256_min_ps(_mm256_max_ps(blended_hi, zero), max_val);

            // 转换回 32-bit 整数
            __m256i result_lo = _mm256_cvtps_epi32(blended_lo);
            __m256i result_hi = _mm256_cvtps_epi32(blended_hi);

            // 打包为 16-bit 整数
            __m256i result_i16 = _mm256_packs_epi32(result_lo, result_hi);
            result_i16 = _mm256_permute4x64_epi64(result_i16, 0xD8);  // 解决乱序问题

            // 打包为 8-bit 无符号整数
            __m256i result_u8 = _mm256_packus_epi16(result_i16, result_i16);
            result_u8 = _mm256_permute4x64_epi64(result_u8, 0xD8);

            // 存储结果（分高低 128-bit）
            if (offset == 0) {
                _mm_storeu_si128(reinterpret_cast<__m128i*>(m_data + i), _mm256_extracti128_si256(result_u8, 0));
            } else {
                _mm_storeu_si128(reinterpret_cast<__m128i*>(m_data + i + 16), _mm256_extracti128_si256(result_u8, 0));
            }
        }
    }

    // 处理剩余像素（标量循环）
    int remainder = total % simd_step;
    int tail_start = total - remainder;
    for (int i = tail_start; i < total; ++i) {
        float blended = alpha * m_data[i] + beta * other.m_data[i];
        m_data[i] = static_cast<unsigned char>(std::clamp(blended, 0.0f, 255.0f));
    }
}

// --- 像素访问 ---
// const unsigned char& Image::at(int x, int y, int channel) const {
//     checkBounds(x, y, channel);
//     return m_data[(y * m_width + x) * m_channels + channel];
// }

} // namespace SimpleImage