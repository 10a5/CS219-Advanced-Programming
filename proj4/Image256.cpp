#include "SimpleImage.h"
#include <algorithm>
#include "stb_image.h"
#include "stb_image_write.h"
#include <string.h>

namespace SimpleImage {

Image256::Image256(int width, int height, int channels) 
    : Image(width, height, channels) {
    if (channels != 1 && channels != 3) {
        throw std::invalid_argument("Image256 only supports 1 or 3 channels");
    }
}

Image256::Image256(const std::string& filename) {
    if (!load(filename)) {
        throw std::runtime_error("Failed to load Image256: " + filename);
    }
}

bool Image256::load(const std::string& filename) {
    int w, h, c;
    unsigned char* data = stbi_load(filename.c_str(), &w, &h, &c, 0);
    if (!data) return false;

    m_width = w;
    m_height = h;
    m_channels = c;
    size_t size = w * h * c;
    m_data = new unsigned char[size];
    memcpy(m_data, data, size);
    stbi_image_free(data);
    
    m_refCount = new int(1);
    return true;
}

bool Image256::save(const std::string& filename) const {
    int success = stbi_write_png(filename.c_str(), 
        static_cast<int>(m_width), 
        static_cast<int>(m_height), 
        m_channels, 
        m_data, 
        m_width * m_channels);
    return success != 0;
}

void Image256::adjustBrightness(int value) {
    value = std::clamp(value, -255, 255);
    const size_t total = m_width * m_height * m_channels;
    for (size_t i = 0; i < total; ++i) {
        int new_val = static_cast<int>(m_data[i]) + value;
        m_data[i] = static_cast<unsigned char>(std::clamp(new_val, 0, 255));
    }
}

void Image256::blend(const Image& other, float alpha) {
    const auto& rhs = dynamic_cast<const Image256&>(other);
    const float beta = 1.0f - alpha;
    const size_t total = m_width * m_height * m_channels;
    for (size_t i = 0; i < total; ++i) {
        float blended = alpha * m_data[i] + beta * rhs.m_data[i];
        m_data[i] = static_cast<unsigned char>(std::clamp(blended, 0.0f, 255.0f));
    }
}

} 