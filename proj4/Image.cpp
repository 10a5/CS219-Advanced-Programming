#include "SimpleImage.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <immintrin.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <bits/chrono.h>
#include <iomanip>

namespace SimpleImage {

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

void Image::checkBounds(int x, int y, int channel) const {
    if (x < 0 || x >= m_width || y < 0 || y >= m_height || channel < 0 || channel >= m_channels) {
        throw std::out_of_range("Pixel access out of bounds");
    }
}

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

Image::~Image() {
    releaseResource();
}

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

bool Image::load(const std::string& filename) {
    releaseResource();
    
    int width, height, channels; 
    unsigned char* stb_data = stbi_load(filename.c_str(), &width, &height, &channels, 0);
    if (!stb_data) {
        throw std::invalid_argument("Failed to load image: " + filename);
        return false;
    }
    if (width < 0 || height < 0) {
        stbi_image_free(stb_data);
        throw std::invalid_argument("Invalid image dimensions (negative width/height)");
    }
    if (static_cast<size_t>(width) > SIZE_MAX / static_cast<size_t>(height)) {
        stbi_image_free(stb_data);
        throw std::invalid_argument("Image dimensions too large");
    }
    m_width = static_cast<size_t>(width);
    m_height = static_cast<size_t>(height);
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
    
    if (m_width > static_cast<size_t>(INT_MAX) || m_height > static_cast<size_t>(INT_MAX)) {
        throw std::invalid_argument("Image dimensions exceed STB's maximum supported size");
    }
    
    const int width = static_cast<int>(m_width);
    const int height = static_cast<int>(m_height);
    
    std::string ext = filename.substr(filename.find_last_of(".") + 1);
    int success = 0;
    
    if (ext == "png") {
        success = stbi_write_png(filename.c_str(), width, height, m_channels, m_data, width * m_channels);
    } else if (ext == "jpg" || ext == "jpeg") {
        success = stbi_write_jpg(filename.c_str(), width, height, m_channels, m_data, 90);
    } else if (ext == "bmp") {
        success = stbi_write_bmp(filename.c_str(), width, height, m_channels, m_data);
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
std::string Image::generate_timestamp_filename(const std::string& prefix) const {
    auto now = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;
    auto timer = std::chrono::system_clock::to_time_t(now);
    
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&timer), "%Y%m%d_%H%M%S")
        << "_" << std::setfill('0') << std::setw(3) << ms.count()
        << "_" << prefix << "." << get_default_extension();
    return oss.str();
}

std::string Image::get_default_extension() const {
    return "png";
}

bool Image::save() const {
    return save(generate_timestamp_filename("image"));
}

void Image::adjustBrightness(int value) {
    if(this->m_refCount == nullptr) {
        throw std::invalid_argument("Uninitiallized image!");
    }
    if (!m_data || value == 0) return;

    const int total = m_width * m_height * m_channels;
    constexpr int simd_step = 32; 

    value = std::clamp(value, -255, 255);
    const __m256i v_value = _mm256_set1_epi16(static_cast<short>(value));

    #pragma omp parallel for
    for (int i = 0; i <= total - simd_step; i += simd_step) {
        __m256i data = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(m_data + i));

        __m256i lo = _mm256_cvtepu8_epi16(_mm256_castsi256_si128(data));  
        __m256i hi = _mm256_cvtepu8_epi16(_mm256_extracti128_si256(data, 1));  
        
        lo = _mm256_adds_epu16(lo, v_value);
        hi = _mm256_adds_epu16(hi, v_value);

        __m256i blended = _mm256_packus_epi16(lo, hi);
        blended = _mm256_permute4x64_epi64(blended, 0xD8);  
        
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(m_data + i), blended);
    }

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
    
    const float beta = 1.0f - alpha;
    const int total = m_width * m_height * m_channels;
    constexpr int simd_step = 32;  
    const __m256 alpha_vec = _mm256_set1_ps(alpha);
    const __m256 beta_vec = _mm256_set1_ps(beta);
    const __m256 zero = _mm256_setzero_ps();
    const __m256 max_val = _mm256_set1_ps(255.0f);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < total - simd_step + 1; i += simd_step) {
        __m256i this_data = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(m_data + i));
        __m256i other_data = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(other.m_data + i));

        for (int offset = 0; offset < 32; offset += 16) {
            __m128i this_chunk = _mm256_extracti128_si256(this_data, offset / 16);
            __m256i this_lo = _mm256_cvtepu8_epi32(this_chunk);
            __m256 this_lo_f = _mm256_cvtepi32_ps(this_lo);
            __m256i this_hi = _mm256_cvtepu8_epi32(_mm_unpackhi_epi64(this_chunk, this_chunk));
            __m256 this_hi_f = _mm256_cvtepi32_ps(this_hi);

            __m128i other_chunk = _mm256_extracti128_si256(other_data, offset / 16);
            __m256i other_lo = _mm256_cvtepu8_epi32(other_chunk);
            __m256 other_lo_f = _mm256_cvtepi32_ps(other_lo);
            __m256i other_hi = _mm256_cvtepu8_epi32(_mm_unpackhi_epi64(other_chunk, other_chunk));
            __m256 other_hi_f = _mm256_cvtepi32_ps(other_hi);

            __m256 blended_lo = _mm256_fmadd_ps(alpha_vec, this_lo_f, _mm256_mul_ps(beta_vec, other_lo_f));
            __m256 blended_hi = _mm256_fmadd_ps(alpha_vec, this_hi_f, _mm256_mul_ps(beta_vec, other_hi_f));

            blended_lo = _mm256_min_ps(_mm256_max_ps(blended_lo, zero), max_val);
            blended_hi = _mm256_min_ps(_mm256_max_ps(blended_hi, zero), max_val);

            __m256i result_lo = _mm256_cvtps_epi32(blended_lo);
            __m256i result_hi = _mm256_cvtps_epi32(blended_hi);

            __m256i result_i16 = _mm256_packs_epi32(result_lo, result_hi);
            result_i16 = _mm256_permute4x64_epi64(result_i16, 0xD8); 

            __m256i result_u8 = _mm256_packus_epi16(result_i16, result_i16);
            result_u8 = _mm256_permute4x64_epi64(result_u8, 0xD8);

            if (offset == 0) {
                _mm_storeu_si128(reinterpret_cast<__m128i*>(m_data + i), _mm256_extracti128_si256(result_u8, 0));
            } else {
                _mm_storeu_si128(reinterpret_cast<__m128i*>(m_data + i + 16), _mm256_extracti128_si256(result_u8, 0));
            }
        }
    }

    int remainder = total % simd_step;
    int tail_start = total - remainder;
    for (int i = tail_start; i < total; ++i) {
        float blended = alpha * m_data[i] + beta * other.m_data[i];
        m_data[i] = static_cast<unsigned char>(std::clamp(blended, 0.0f, 255.0f));
    }
}

} 