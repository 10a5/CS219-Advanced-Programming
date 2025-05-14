#include "SimpleImage.h"
#include <algorithm>
#include "stb_image_write.h"

namespace SimpleImage {

JpegImage::JpegImage(int width, int height) 
    : Image(width, height, 3) {} 

JpegImage::JpegImage(const std::string& filename) {
    if (!load(filename)) {
        throw std::runtime_error("Failed to load JPEG: " + filename);
    }
}

bool JpegImage::save(const std::string& filename) const {
    if (!m_data) return false;
    int width = static_cast<int>(m_width);
    int height = static_cast<int>(m_height);
    return stbi_write_jpg(filename.c_str(), width, height, m_channels, m_data, m_quality);
}

void JpegImage::setQuality(int quality) {
    m_quality = std::clamp(quality, 1, 100);
}
} 