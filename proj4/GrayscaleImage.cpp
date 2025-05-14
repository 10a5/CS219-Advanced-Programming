#include "SimpleImage.h"
#include <algorithm>

namespace SimpleImage {

GrayscaleImage::GrayscaleImage(int width, int height) 
    : Image(width, height, 1) {} 

GrayscaleImage::GrayscaleImage(const std::string& filename) {
    if (!load(filename)) {
        throw std::runtime_error("Failed to load grayscale image: " + filename);
    }
}

bool GrayscaleImage::load(const std::string& filename) {
    bool success = Image::load(filename);
    if (success && m_channels > 1) {
        unsigned char* grayData = new unsigned char[m_width * m_height];
        for (size_t i = 0; i < m_width * m_height * m_channels; i += m_channels) {
            grayData[i / m_channels] = static_cast<unsigned char>(
                0.299f * m_data[i] + 0.587f * m_data[i+1] + 0.114f * m_data[i+2]);
        }
        delete[] m_data;
        m_data = grayData;
        m_channels = 1;
    }
    return success;
}

unsigned char GrayscaleImage::getPixel(int x, int y) const {
    checkBounds(x, y, 0);
    return m_data[y * m_width + x];
}

void GrayscaleImage::setPixel(int x, int y, unsigned char value) {
    checkBounds(x, y, 0);
    m_data[y * m_width + x] = value;
}
} 