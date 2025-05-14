#include "SimpleImage.h"
#include "stb_image_write.h"

namespace SimpleImage {

PngImage::PngImage(int width, int height, int channels) 
    : Image(width, height, channels) {
    if (channels != 3 && channels != 4) {
        throw std::invalid_argument("PNG must have 3 or 4 channels");
    }
}

PngImage::PngImage(const std::string& filename) {
    if (!load(filename)) {
        throw std::runtime_error("Failed to load PNG: " + filename);
    }
}
}