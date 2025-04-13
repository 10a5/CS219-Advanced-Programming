#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <immintrin.h>
#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#pragma pack(push, 1)
typedef struct {
    uint16_t bfType;
    uint32_t bfSize;
    uint16_t bfReserved1;
    uint16_t bfReserved2;
    uint32_t bfOffBits;
} BITMAPFILEHEADER;

typedef struct {
    uint32_t biSize;
    int32_t biWidth;
    int32_t biHeight;
    uint16_t biPlanes;
    uint16_t biBitCount;
    uint32_t biCompression;
    uint32_t biSizeImage;
    int32_t biXPelsPerMeter;
    int32_t biYPelsPerMeter;
    uint32_t biClrUsed;
    uint32_t biClrImportant;
} BITMAPINFOHEADER;
#pragma pack(pop)

typedef struct {
    BITMAPFILEHEADER file_header;
    BITMAPINFOHEADER info_header;
    int width;
    int height;
    unsigned char *data;
} BMPImage;

void free_bmp(BMPImage *image) {
    if (image) {
        free(image->data);
        free(image);
    }
}

BMPImage *read_bmp(const char *filename, char **error) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        *error = "Failed to open file";
        return NULL;
    }

    BITMAPFILEHEADER file_header;
    if (fread(&file_header, sizeof(file_header), 1, file) != 1) {
        *error = "Failed to read file header";
        fclose(file);
        return NULL;
    }

    if (file_header.bfType != 0x4D42) {
        *error = "Not a BMP file";
        fclose(file);
        return NULL;
    }

    BITMAPINFOHEADER info_header;
    if (fread(&info_header, sizeof(info_header), 1, file) != 1) {
        *error = "Failed to read info header";
        fclose(file);
        return NULL;
    }

    if (info_header.biBitCount != 24 || info_header.biCompression != 0) {
        *error = "Only 24-bit uncompressed BMP supported";
        fclose(file);
        return NULL;
    }

    const int src_height = abs(info_header.biHeight);
    const int width = info_header.biWidth;
    const size_t row_size = ((width * 3 + 3) / 4) * 4;
    const size_t data_size = row_size * src_height;

    unsigned char *pixel_data = malloc(data_size);
    if (!pixel_data) {
        *error = "Memory allocation failed";
        fclose(file);
        return NULL;
    }

    fseek(file, file_header.bfOffBits, SEEK_SET);
    if (fread(pixel_data, 1, data_size, file) != data_size) {
        *error = "Failed to read pixel data";
        free(pixel_data);
        fclose(file);
        return NULL;
    }
    fclose(file);

    const size_t plane_size = width * src_height;
    unsigned char *planar_data = malloc(plane_size * 3);
    if (!planar_data) {
        *error = "Memory allocation failed";
        free(pixel_data);
        return NULL;
    }

    unsigned char *r_plane = planar_data;
    unsigned char *g_plane = planar_data + plane_size;
    unsigned char *b_plane = planar_data + plane_size * 2;

    if (info_header.biHeight > 0) {
        for (int y = 0; y < src_height; y++) {
            const unsigned char *src_row = pixel_data + (src_height - 1 - y) * row_size;
            for (int x = 0; x < width; x++) {
                const int dst_idx = y * width + x;
                b_plane[dst_idx] = src_row[x * 3];
                g_plane[dst_idx] = src_row[x * 3 + 1];
                r_plane[dst_idx] = src_row[x * 3 + 2];
            }
        }
    } else {
        for (int y = 0; y < src_height; y++) {
            const unsigned char *src_row = pixel_data + y * row_size;
            for (int x = 0; x < width; x++) {
                const int dst_idx = y * width + x;
                b_plane[dst_idx] = src_row[x * 3];
                g_plane[dst_idx] = src_row[x * 3 + 1];
                r_plane[dst_idx] = src_row[x * 3 + 2];
            }
        }
    }
    free(pixel_data);

    BMPImage *image = malloc(sizeof(BMPImage));
    image->file_header = file_header;
    image->info_header = info_header;
    image->width = width;
    image->height = src_height;
    image->data = planar_data;
    image->info_header.biHeight = -src_height;
    return image;
}

int write_bmp(const char *filename, BMPImage *image, char **error) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        *error = "Failed to create file";
        return 0;
    }

    const size_t row_size = ((image->width * 3 + 3) / 4) * 4;
    const size_t data_size = row_size * image->height;
    unsigned char *padded_data = calloc(1, data_size);
    if (!padded_data) {
        *error = "Memory allocation failed";
        fclose(file);
        return 0;
    }

    const size_t plane_size = image->width * image->height;
    const unsigned char *r_plane = image->data;
    const unsigned char *g_plane = image->data + plane_size;
    const unsigned char *b_plane = image->data + plane_size * 2;

    for (int y = 0; y < image->height; y++) {
        unsigned char *dst_row = padded_data + y * row_size;
        for (int x = 0; x < image->width; x++) {
            const int src_idx = y * image->width + x;
            dst_row[x * 3] = b_plane[src_idx];
            dst_row[x * 3 + 1] = g_plane[src_idx];
            dst_row[x * 3 + 2] = r_plane[src_idx];
        }
    }

    image->file_header.bfSize = sizeof(BITMAPFILEHEADER) + 
                               sizeof(BITMAPINFOHEADER) + 
                               data_size;
    image->info_header.biSizeImage = data_size;
    image->file_header.bfOffBits = sizeof(BITMAPFILEHEADER) + 
                                  sizeof(BITMAPINFOHEADER);
    image->info_header.biHeight = -image->height;

    int success = 1;
    if (fwrite(&image->file_header, sizeof(BITMAPFILEHEADER), 1, file) != 1 ||
        fwrite(&image->info_header, sizeof(BITMAPINFOHEADER), 1, file) != 1 ||
        fwrite(padded_data, 1, data_size, file) != data_size) {
        *error = "Write failed";
        success = 0;
    }

    free(padded_data);
    fclose(file);
    return success;
}
void adjust_brightness(BMPImage *image, int value) {
    const size_t total = 3 * image->width * image->height;
    unsigned char *data = image->data;
    const int vector_size = 16;
    const size_t end = total - vector_size;

    #pragma omp parallel
    {
        #pragma omp for
        for (size_t i = 0; i <= end; i += vector_size) {
            __m128i vec = _mm_loadu_si128((__m128i*)(data + i));
            
            __m128i low = _mm_cvtepu8_epi16(vec);
            __m128i high = _mm_cvtepu8_epi16(_mm_srli_si128(vec, 8));
            
            __m128i val_vec = _mm_set1_epi16(value);
            
            low = _mm_add_epi16(low, val_vec);
            high = _mm_add_epi16(high, val_vec);
        
            low = _mm_min_epi16(_mm_max_epi16(low, _mm_setzero_si128()), _mm_set1_epi16(255));
            high = _mm_min_epi16(_mm_max_epi16(high, _mm_setzero_si128()), _mm_set1_epi16(255));
            
            __m128i res = _mm_packus_epi16(low, high);
            
            _mm_storeu_si128((__m128i*)(data + i), res);
        }

        #pragma omp for
        for (size_t i = ((total / vector_size) * vector_size); i < total; ++i) {
            int temp = data[i] + value;
            data[i] = (unsigned char)(temp > 255 ? 255 : (temp < 0 ? 0 : temp));
        }
    }
}

BMPImage *blend_images(BMPImage *img1, BMPImage *img2, char **error) {
    if (img1->width != img2->width || img1->height != img2->height) {
        *error = "Graphic size mismatch";
        return NULL;
    }

    BMPImage *result = malloc(sizeof(BMPImage));
    if (!result) {
        *error = "Fail to allign memory";
        return NULL;
    }

    *result = *img1;
    const int total_pixels = img1->width * img1->height;
    const int total_bytes = total_pixels * 3;
    
    result->data = malloc(total_bytes);
    if (!result->data) {
        *error = "Fail to allign memory";
        free(result);
        return NULL;
    }

    unsigned char *p1 = img1->data;
    unsigned char *p2 = img2->data;
    unsigned char *pr = result->data;

    #pragma omp parallel for simd schedule(simd:static) aligned(pr, p1, p2:32)
    for (int i = 0; i < total_bytes; i += 32) {
        __m256i data1 = _mm256_loadu_si256((__m256i*)(p1 + i));
        __m256i data2 = _mm256_loadu_si256((__m256i*)(p2 + i));

        __m256i avg = _mm256_avg_epu8(data1, data2);

        _mm256_storeu_si256((__m256i*)(pr + i), avg);
    }

    #pragma omp parallel for
    for (int i = (total_bytes / 32) * 32; i < total_bytes; ++i) {
        pr[i] = (p1[i] + p2[i]) / 2;
    }

    return result;
}

void adjust_contrast(BMPImage *image, float alpha) {
    unsigned char *p = image->data;
    const int total_bytes = image->width * image->height * 3;
    
    const int alpha_q15 = (int)(alpha * 16384);
    const __m256i _128_vec = _mm256_set1_epi16(128);

    const int aligned_bytes = total_bytes - (total_bytes % 32);

    const __m256i alpha_vec = _mm256_set1_epi16(alpha_q15);

    #pragma omp parallel for
    for (int i = 0; i < aligned_bytes; i += 32) {
        __m256i data = _mm256_loadu_si256((__m256i*)(p + i));
        __m256i lo = _mm256_cvtepu8_epi16(_mm256_extracti128_si256(data, 0));
        __m256i hi = _mm256_cvtepu8_epi16(_mm256_extracti128_si256(data, 1));

        lo = _mm256_sub_epi16(lo, _128_vec);
        hi = _mm256_sub_epi16(hi, _128_vec);

        lo = _mm256_slli_epi16(lo, 8);
        hi = _mm256_slli_epi16(hi, 8);

        lo = _mm256_mulhrs_epi16(lo, alpha_vec);
        hi = _mm256_mulhrs_epi16(hi, alpha_vec);

        lo = _mm256_srai_epi16(lo, 8); 
        hi = _mm256_srai_epi16(hi, 8);

        lo = _mm256_add_epi16(lo, lo);
        hi = _mm256_add_epi16(hi, hi);
        lo = _mm256_add_epi16(lo, _128_vec);
        hi = _mm256_add_epi16(hi, _128_vec);

        __m256i packed = _mm256_packus_epi16(lo, hi);
        packed = _mm256_permute4x64_epi64(packed, _MM_SHUFFLE(3,1,2,0));
        _mm256_storeu_si256((__m256i*)(p + i), packed);
    }

    #pragma omp parallel for
    for (int i = aligned_bytes; i < total_bytes; ++i) {
        int val = ((p[i] - 128) * alpha_q15 + 0x4000) >> 15; 
        val += 128;
        p[i] = (unsigned char)(val > 255 ? 255 : (val < 0 ? 0 : val));
    }
}

void convert_to_grayscale(BMPImage *image) {
    const size_t plane_size = image->width * image->height;
    const size_t total_pixels = plane_size;
    unsigned char *r = image->data;
    unsigned char *g = r + plane_size;
    unsigned char *b = g + plane_size;

    const __m256 r_coeff = _mm256_set1_ps(0.299f);
    const __m256 g_coeff = _mm256_set1_ps(0.587f);
    const __m256 b_coeff = _mm256_set1_ps(0.114f);
    const __m256 max_val = _mm256_set1_ps(255.0f);
    const __m256 zero = _mm256_setzero_ps();

    #pragma omp parallel for schedule(guided)
    for (size_t i = 0; i < total_pixels; i += 8) {
        __m128i r8 = _mm_loadl_epi64((__m128i*)(r + i));
        __m128i g8 = _mm_loadl_epi64((__m128i*)(g + i));
        __m128i b8 = _mm_loadl_epi64((__m128i*)(b + i));

        __m256i r32 = _mm256_cvtepu8_epi32(r8);
        __m256i g32 = _mm256_cvtepu8_epi32(g8);
        __m256i b32 = _mm256_cvtepu8_epi32(b8);

        __m256 rf = _mm256_cvtepi32_ps(r32);
        __m256 gf = _mm256_cvtepi32_ps(g32);
        __m256 bf = _mm256_cvtepi32_ps(b32);

        __m256 gray_f = _mm256_add_ps(
            _mm256_mul_ps(rf, r_coeff),
            _mm256_mul_ps(gf, g_coeff));
        gray_f = _mm256_add_ps(gray_f, _mm256_mul_ps(bf, b_coeff));

        gray_f = _mm256_min_ps(_mm256_max_ps(gray_f, zero), max_val);

        __m256i gray32 = _mm256_cvtps_epi32(gray_f);
        __m128i gray16 = _mm_packus_epi32(
            _mm256_extractf128_si256(gray32, 0),
            _mm256_extractf128_si256(gray32, 1));
        __m128i gray8 = _mm_packus_epi16(gray16, gray16);

        _mm_storel_epi64((__m128i*)(r + i), gray8);
        _mm_storel_epi64((__m128i*)(g + i), gray8);
        _mm_storel_epi64((__m128i*)(b + i), gray8);
    }

    #pragma omp parallel for
    for (size_t i = (total_pixels / 8) * 8; i < total_pixels; ++i) {
        const float gray = 
            r[i] * 0.299f + 
            g[i] * 0.587f + 
            b[i] * 0.114f;
        
        const unsigned char gray_val = (unsigned char)fminf(fmaxf(gray, 0.0f), 255.0f);
        r[i] = g[i] = b[i] = gray_val;
    }
}

void invert_colors(BMPImage *image) {
    const size_t total = 3 * image->width * image->height;
    unsigned char *data = image->data;
    const __m128i ones = _mm_set1_epi8(0xFF);
    const int vector_size = 16;  
    const size_t end = total - vector_size;

    #pragma omp parallel
    {
        #pragma omp for
        for (size_t i = 0; i <= end; i += vector_size) {
            __m128i vec = _mm_loadu_si128((__m128i*)(data + i));
            
            __m128i res = _mm_xor_si128(vec, ones);
            
            _mm_storeu_si128((__m128i*)(data + i), res);
        }

        #pragma omp for
        for (size_t i = ((total / vector_size) * vector_size); i < total; ++i) {
            data[i] = 255 - data[i];
        }
    }
}

BMPImage *crop_image(BMPImage *image, int x, int y, int w, int h, char **error) {
    if (x < 0 || y < 0 || x + w > image->width || y + h > image->height) {
        *error = "Crop area exceeds image bounds";
        return NULL;
    }

    BMPImage *result = malloc(sizeof(BMPImage));
    *result = *image;
    result->width = w;
    result->height = h;
    result->info_header.biWidth = w;
    result->info_header.biHeight = h;

    size_t new_plane_size = w * h;
    result->data = malloc(new_plane_size * 3);
    if (!result->data) {
        *error = "Failed to allocate memory";
        free(result);
        return NULL;
    }

    unsigned char *src_r = image->data;
    unsigned char *src_g = src_r + image->width * image->height;
    unsigned char *src_b = src_g + image->width * image->height;

    unsigned char *dst_r = result->data;
    unsigned char *dst_g = dst_r + new_plane_size;
    unsigned char *dst_b = dst_g + new_plane_size;

    #pragma omp parallel for schedule(guided)
    for (int row = 0; row < h; row++) {
        const int src_y = y + row;
        if (src_y >= 0 && src_y < image->height) {
            for (int col = 0; col < w; col++) {
                const int src_x = x + col;
                if (src_x >= 0 && src_x < image->width) {
                    const int src_idx = src_y * image->width + src_x;
                    const int dst_idx = row * w + col;
                    
                    dst_r[dst_idx] = src_r[src_idx];
                    dst_g[dst_idx] = src_g[src_idx];
                    dst_b[dst_idx] = src_b[src_idx];
                }
            }
        }
    }
    return result;
}
BMPImage *rotate_image(BMPImage *image, int degrees, char **error) {
    if (degrees != 90 && degrees != -90 && degrees != 180) {
        *error = "Only 90, -90, and 180 degree rotations are supported";
        return NULL;
    }

    BMPImage *result = malloc(sizeof(BMPImage));
    *result = *image;

    const int is_90_degree = (degrees == 90 || degrees == -90);
    if (is_90_degree) {
        result->width = image->height;
        result->height = image->width;
        result->info_header.biWidth = image->height;
        result->info_header.biHeight = image->width;
    }

    const size_t new_plane_size = result->width * result->height;
    unsigned char *planar_data = aligned_alloc(64, new_plane_size * 3);
    if (!planar_data) {
        *error = "Failed to allocate memory";
        free(result);
        return NULL;
    }
    result->data = planar_data;

    const unsigned char *src_r = image->data;
    const unsigned char *src_g = src_r + image->width * image->height;
    const unsigned char *src_b = src_g + image->width * image->height;

    unsigned char *dst_r = planar_data;
    unsigned char *dst_g = dst_r + new_plane_size;
    unsigned char *dst_b = dst_g + new_plane_size;

    switch (degrees) {
        case 90:
            #pragma omp parallel for schedule(static) collapse(2)
            for (int y = 0; y < image->height; y++) {
                for (int x = 0; x < image->width; x++) {
                    const int dst_x = image->height - 1 - y;
                    const int dst_y = x;
                    const int dst_idx = dst_y * result->width + dst_x;
                    const int src_idx = y * image->width + x;
                    
                    dst_r[dst_idx] = src_r[src_idx];
                    dst_g[dst_idx] = src_g[src_idx];
                    dst_b[dst_idx] = src_b[src_idx];
                }
            }
            break;

        case -90:
            #pragma omp parallel for schedule(static) collapse(2)
            for (int y = 0; y < image->height; y++) {
                for (int x = 0; x < image->width; x++) {
                    const int dst_x = y;
                    const int dst_y = image->width - 1 - x;
                    const int dst_idx = dst_y * result->width + dst_x;
                    const int src_idx = y * image->width + x;
                    
                    dst_r[dst_idx] = src_r[src_idx];
                    dst_g[dst_idx] = src_g[src_idx];
                    dst_b[dst_idx] = src_b[src_idx];
                }
            }
            break;

        case 180:
            #pragma omp parallel for schedule(static)
            for (int y = 0; y < image->height; y++) {
                for (int x = 0; x < image->width; x++) {
                    const int dst_x = image->width - 1 - x;
                    const int dst_y = image->height - 1 - y;
                    const int dst_idx = dst_y * result->width + dst_x;
                    const int src_idx = y * image->width + x;
                    
                    dst_r[dst_idx] = src_r[src_idx];
                    dst_g[dst_idx] = src_g[src_idx];
                    dst_b[dst_idx] = src_b[src_idx];
                }
            }
            break;
    }

    return result;
}

BMPImage *resize_image(BMPImage *image, int new_width, int new_height, char **error) {
    if (new_width <= 0 || new_height <= 0) {
        *error = "The resize dimensions must be greater than zero";
        return NULL;
    }

    BMPImage *result = malloc(sizeof(BMPImage));
    *result = *image;
    result->width = new_width;
    result->height = new_height;
    result->info_header.biWidth = new_width;
    result->info_header.biHeight = new_height;

    const size_t new_plane_size = new_width * new_height;
    unsigned char *planar_data = aligned_alloc(64, new_plane_size * 3);
    if (!planar_data) {
        *error = "Failed to allocate memory";
        free(result);
        return NULL;
    }
    result->data = planar_data;

    const unsigned char *src_r = image->data;
    const unsigned char *src_g = src_r + image->width * image->height;
    const unsigned char *src_b = src_g + image->width * image->height;

    unsigned char *dst_r = planar_data;
    unsigned char *dst_g = dst_r + new_plane_size;
    unsigned char *dst_b = dst_g + new_plane_size;

    const float x_ratio = (float)image->width / new_width;
    const float y_ratio = (float)image->height / new_height;

    #pragma omp parallel for schedule(guided) collapse(2)
    for (int y = 0; y < new_height; y++) {
        for (int x = 0; x < new_width; x++) {
            const int src_x = (int)(x * x_ratio);
            const int src_y = (int)(y * y_ratio);
            const int src_idx = (src_y < image->height ? src_y : image->height-1) * image->width 
                              + (src_x < image->width ? src_x : image->width-1);
            const int dst_idx = y * new_width + x;

            dst_r[dst_idx] = src_r[src_idx];
            dst_g[dst_idx] = src_g[src_idx];
            dst_b[dst_idx] = src_b[src_idx];
        }
    }

    return result;
}

void apply_vignette(BMPImage *image, float intensity) {
    const int width = image->width;
    const int height = image->height;
    const size_t plane_size = width * height;
    unsigned char *r = image->data;
    unsigned char *g = image->data + plane_size;
    unsigned char *b = image->data + plane_size * 2;

    const float center_x = (width - 1) * 0.5f;
    const float center_y = (height - 1) * 0.5f;
    const float max_radius_sq = center_x*center_x + center_y*center_y;
    const __m256 intensity_vec = _mm256_set1_ps(intensity / max_radius_sq);
    const __m256 one_vec = _mm256_set1_ps(1.0f);
    const __m256 zero_vec = _mm256_setzero_ps();
    const __m256i mask = _mm256_setr_epi32(0,1,2,3,4,5,6,7);

    #pragma omp parallel for schedule(guided)
    for (int y = 0; y < height; ++y) {
        const __m256 dy = _mm256_set1_ps(y - center_y);
        const __m256 dy_sq = _mm256_mul_ps(dy, dy);
        unsigned char *r_row = r + y*width;
        unsigned char *g_row = g + y*width;
        unsigned char *b_row = b + y*width;

        int x = 0;
        for (; x <= width-8; x += 8) {
            __m256i x_idx = _mm256_add_epi32(_mm256_set1_epi32(x), mask);
            __m256 x_f = _mm256_cvtepi32_ps(x_idx);
            
            __m256 dx = _mm256_sub_ps(x_f, _mm256_set1_ps(center_x));
            __m256 dx_sq = _mm256_mul_ps(dx, dx);
            
            __m256 distance_sq = _mm256_add_ps(dx_sq, dy_sq);
            
            __m256 factor = _mm256_mul_ps(distance_sq, intensity_vec);
            factor = _mm256_sub_ps(one_vec, factor);

            __m256i r_data = _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i*)(r_row+x)));
            __m256i g_data = _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i*)(g_row+x)));
            __m256i b_data = _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i*)(b_row+x)));

            __m256 r_f = _mm256_cvtepi32_ps(r_data);
            __m256 g_f = _mm256_cvtepi32_ps(g_data);
            __m256 b_f = _mm256_cvtepi32_ps(b_data);

            r_f = _mm256_mul_ps(r_f, factor);
            g_f = _mm256_mul_ps(g_f, factor);
            b_f = _mm256_mul_ps(b_f, factor);

            r_f = _mm256_max_ps(_mm256_min_ps(r_f, _mm256_set1_ps(255.0f)), zero_vec);
            g_f = _mm256_max_ps(_mm256_min_ps(g_f, _mm256_set1_ps(255.0f)), zero_vec);
            b_f = _mm256_max_ps(_mm256_min_ps(b_f, _mm256_set1_ps(255.0f)), zero_vec);

            __m256i r_result = _mm256_cvtps_epi32(r_f);
            __m256i g_result = _mm256_cvtps_epi32(g_f);
            __m256i b_result = _mm256_cvtps_epi32(b_f);

            __m128i low_r = _mm256_castsi256_si128(r_result);
            __m128i high_r = _mm256_extracti128_si256(r_result, 1);
            __m128i low_g = _mm256_castsi256_si128(g_result);
            __m128i high_g = _mm256_extracti128_si256(g_result, 1);
            __m128i low_b = _mm256_castsi256_si128(b_result);
            __m128i high_b = _mm256_extracti128_si256(b_result, 1);

            const __m128i shuffle_mask = _mm_setr_epi8(
                0, 1, 4, 5, 8, 9, 12, 13, 
                0, 1, 4, 5, 8, 9, 12, 13 
            );

            __m128i packed_low_r = _mm_shuffle_epi8(low_r, shuffle_mask); 
            __m128i packed_high_r = _mm_shuffle_epi8(high_r, shuffle_mask);
            __m128i packed_low_g = _mm_shuffle_epi8(low_g, shuffle_mask); 
            __m128i packed_high_g = _mm_shuffle_epi8(high_g, shuffle_mask);
            __m128i packed_low_b = _mm_shuffle_epi8(low_b, shuffle_mask);
            __m128i packed_high_b = _mm_shuffle_epi8(high_b, shuffle_mask);
            __m128i result_r = _mm_unpacklo_epi64(packed_low_r, packed_high_r);
            __m128i result_g = _mm_unpacklo_epi64(packed_low_g, packed_high_g);
            __m128i result_b = _mm_unpacklo_epi64(packed_low_b, packed_high_b);

            const __m128i shuffle_mask1 = _mm_setr_epi8(
                0, 2, 4, 6, 8, 10, 12, 14,
                -1, -1, -1, -1, -1, -1, -1, -1 
            );
            __m128i r_r = _mm_shuffle_epi8(result_r, shuffle_mask1);
            __m128i g_r = _mm_shuffle_epi8(result_g, shuffle_mask1);
            __m128i b_r = _mm_shuffle_epi8(result_b, shuffle_mask1);

            _mm_storel_epi64((__m128i*)(r_row+x), r_r);
            _mm_storel_epi64((__m128i*)(g_row+x), g_r);
            _mm_storel_epi64((__m128i*)(b_row+x), b_r);
        }

        for (; x < width; ++x) {
            const float dx = x - center_x;
            const float distance_sq = dx*dx + (y-center_y)*(y-center_y);
            const float factor = 1.0f - intensity * distance_sq / max_radius_sq;
            
            const unsigned char r_val = r_row[x] * factor;
            const unsigned char g_val = g_row[x] * factor;
            const unsigned char b_val = b_row[x] * factor;
            
            r_row[x] = r_val > 255 ? 255 : (r_val < 0 ? 0 : r_val);
            g_row[x] = g_val > 255 ? 255 : (g_val < 0 ? 0 : g_val);
            b_row[x] = b_val > 255 ? 255 : (b_val < 0 ? 0 : b_val);
        }
    }
}

void adjust_rgb_channels(BMPImage *image, float r_scale, float g_scale, float b_scale) {
    const size_t plane_size = image->width * image->height;
    const __m256 r_scale_vec = _mm256_set1_ps(r_scale);
    const __m256 g_scale_vec = _mm256_set1_ps(g_scale);
    const __m256 b_scale_vec = _mm256_set1_ps(b_scale);
    const __m256 max_val = _mm256_set1_ps(255.0f);
    const __m256 min_val = _mm256_setzero_ps();

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            unsigned char *r = image->data;
            #pragma omp parallel for simd
            for (size_t i = 0; i < plane_size; i += 8) {
                __m128i vec8 = _mm_loadl_epi64((__m128i*)(r + i));
                __m256i vec32 = _mm256_cvtepu8_epi32(vec8);
                __m256 fvec = _mm256_cvtepi32_ps(vec32);
                fvec = _mm256_mul_ps(fvec, r_scale_vec);
                fvec = _mm256_min_ps(_mm256_max_ps(fvec, min_val), max_val);
                __m256i res32 = _mm256_cvtps_epi32(fvec);
                __m128i res16 = _mm_packus_epi32(_mm256_extracti128_si256(res32, 0),
                                               _mm256_extracti128_si256(res32, 1));
                __m128i res8 = _mm_packus_epi16(res16, res16);
                _mm_storel_epi64((__m128i*)(r + i), res8);
            }
        }

        #pragma omp section
        {
            unsigned char *g = image->data + plane_size;
            #pragma omp parallel for simd
            for (size_t i = 0; i < plane_size; i += 8) {
                __m128i vec8 = _mm_loadl_epi64((__m128i*)(g + i));
                __m256i vec32 = _mm256_cvtepu8_epi32(vec8);
                __m256 fvec = _mm256_cvtepi32_ps(vec32);
                fvec = _mm256_mul_ps(fvec, g_scale_vec);  
                fvec = _mm256_min_ps(_mm256_max_ps(fvec, min_val), max_val);
                __m256i res32 = _mm256_cvtps_epi32(fvec);
                __m128i res16 = _mm_packus_epi32(_mm256_extracti128_si256(res32, 0),
                                               _mm256_extracti128_si256(res32, 1));
                __m128i res8 = _mm_packus_epi16(res16, res16);
                _mm_storel_epi64((__m128i*)(g + i), res8);
            }
        }

        #pragma omp section
        {
            unsigned char *b = image->data + plane_size * 2;
            #pragma omp parallel for simd
            for (size_t i = 0; i < plane_size; i += 8) {
                __m128i vec8 = _mm_loadl_epi64((__m128i*)(b + i));
                __m256i vec32 = _mm256_cvtepu8_epi32(vec8);
                __m256 fvec = _mm256_cvtepi32_ps(vec32);
                fvec = _mm256_mul_ps(fvec, b_scale_vec); 
                fvec = _mm256_min_ps(_mm256_max_ps(fvec, min_val), max_val);
                __m256i res32 = _mm256_cvtps_epi32(fvec);
                __m128i res16 = _mm_packus_epi32(_mm256_extracti128_si256(res32, 0),
                                               _mm256_extracti128_si256(res32, 1));
                __m128i res8 = _mm_packus_epi16(res16, res16);
                _mm_storel_epi64((__m128i*)(b + i), res8);
            }
        }
    }

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            unsigned char *r = image->data;
            for (size_t i = (plane_size / 8) * 8; i < plane_size; ++i) {
                float val = r[i] * r_scale;
                r[i] = (unsigned char)fminf(fmaxf(val, 0.0f), 255.0f);
            }
        }
        
        #pragma omp section
        {
            unsigned char *g = image->data + plane_size;
            for (size_t i = (plane_size / 8) * 8; i < plane_size; ++i) {
                float val = g[i] * g_scale;
                g[i] = (unsigned char)fminf(fmaxf(val, 0.0f), 255.0f);
            }
        }
        
        #pragma omp section
        {
            unsigned char *b = image->data + plane_size * 2;
            for (size_t i = (plane_size / 8) * 8; i < plane_size; ++i) {
                float val = b[i] * b_scale;
                b[i] = (unsigned char)fminf(fmaxf(val, 0.0f), 255.0f);
            }
        }
    }
}

void partial_blur(BMPImage *image, int x, int y, int w, int h, int radius) {
    const int width = image->width;
    const int height = image->height;
    const size_t plane_size = width * height;

    unsigned char* temp = aligned_alloc(64, plane_size * 3);
    memcpy(temp, image->data, plane_size * 3);
    unsigned char* temp_r = temp;
    unsigned char* temp_g = temp + plane_size;
    unsigned char* temp_b = temp + plane_size * 2;

    unsigned char* r = image->data;
    unsigned char* g = r + plane_size;
    unsigned char* b = g + plane_size;

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            #pragma omp parallel for collapse(2) schedule(dynamic)
            for (int j = y; j < y + h; j++) {
                for (int i = x; i < x + w; i++) {
                    int total = 0;
                    int count = 0;
                    
                    for (int dy = -radius; dy <= radius; dy++) {
                        const int ny = j + dy;
                        if (ny < 0 || ny >= height) { 
                            continue;
                        }
                        
                        const int start_x = MAX(i - radius, 0);
                        const int end_x = MIN(i + radius, width - 1);
                        const int len = end_x - start_x + 1;
                        
                        __m256i vec_sum = _mm256_setzero_si256();
                        int k = 0;
                        for (; k <= len - 8; k += 8) {
                            const __m128i pixels = _mm_loadl_epi64(
                                (__m128i*)(temp_r + ny * width + start_x + k));
                            const __m256i pixels32 = _mm256_cvtepu8_epi32(pixels);
                            vec_sum = _mm256_add_epi32(vec_sum, pixels32);
                        }
                        
                        __m128i sum128 = _mm_add_epi32(
                            _mm256_extracti128_si256(vec_sum, 0),
                            _mm256_extracti128_si256(vec_sum, 1));
                        sum128 = _mm_add_epi32(sum128, _mm_shuffle_epi32(sum128, _MM_SHUFFLE(1,0,3,2)));
                        sum128 = _mm_add_epi32(sum128, _mm_shuffle_epi32(sum128, _MM_SHUFFLE(2,3,0,1)));
                        total += _mm_extract_epi32(sum128, 0);

                        for (; k < len; k++) {
                            total += temp_r[ny * width + start_x + k];
                        }
                        count += len;
                    }
                
                    if (count > 0) {
                        r[j * width + i] = (total + count/2) / count; 
                    }
                }
            }
        }

        #pragma omp section
        {
            #pragma omp parallel for collapse(2) schedule(dynamic)
            for (int j = y; j < y + h; j++) {
                for (int i = x; i < x + w; i++) {
                    int total = 0;
                    int count = 0;
                    
                    for (int dy = -radius; dy <= radius; dy++) {
                        const int ny = j + dy;
                        if (ny < 0 || ny >= height) {
                            continue;
                        }
                        
                        const int start_x = MAX(i - radius, 0);
                        const int end_x = MIN(i + radius, width - 1);
                        const int len = end_x - start_x + 1;
                        
                        __m256i vec_sum = _mm256_setzero_si256();
                        int k = 0;
                        for (; k <= len - 8; k += 8) {
                            const __m128i pixels = _mm_loadl_epi64(
                                (__m128i*)(temp_g + ny * width + start_x + k));
                            const __m256i pixels32 = _mm256_cvtepu8_epi32(pixels);
                            vec_sum = _mm256_add_epi32(vec_sum, pixels32);
                        }
                        
                        __m128i sum128 = _mm_add_epi32(
                            _mm256_extracti128_si256(vec_sum, 0),
                            _mm256_extracti128_si256(vec_sum, 1));
                        sum128 = _mm_add_epi32(sum128, _mm_shuffle_epi32(sum128, _MM_SHUFFLE(1,0,3,2)));
                        sum128 = _mm_add_epi32(sum128, _mm_shuffle_epi32(sum128, _MM_SHUFFLE(2,3,0,1)));
                        total += _mm_extract_epi32(sum128, 0);

                        for (; k < len; k++) {
                            total += temp_g[ny * width + start_x + k];
                        }
                        count += len;
                    }
                    
                    if (count > 0) {
                        g[j * width + i] = (total + count/2) / count;
                    }
                }
            }
        }

        #pragma omp section
        {
            #pragma omp parallel for collapse(2) schedule(dynamic)
            for (int j = y; j < y + h; j++) {
                for (int i = x; i < x + w; i++) {
                    int total = 0;
                    int count = 0;
                    
                    for (int dy = -radius; dy <= radius; dy++) {
                        const int ny = j + dy;
                        if (ny < 0 || ny >= height) {
                            continue;
                        }
                        
                        const int start_x = MAX(i - radius, 0);
                        const int end_x = MIN(i + radius, width - 1);
                        const int len = end_x - start_x + 1;
                        
                        __m256i vec_sum = _mm256_setzero_si256();
                        int k = 0;
                        for (; k <= len - 8; k += 8) {
                            const __m128i pixels = _mm_loadl_epi64(
                                (__m128i*)(temp_b + ny * width + start_x + k));
                            const __m256i pixels32 = _mm256_cvtepu8_epi32(pixels);
                            vec_sum = _mm256_add_epi32(vec_sum, pixels32);
                        }
                        
                        __m128i sum128 = _mm_add_epi32(
                            _mm256_extracti128_si256(vec_sum, 0),
                            _mm256_extracti128_si256(vec_sum, 1));
                        sum128 = _mm_add_epi32(sum128, _mm_shuffle_epi32(sum128, _MM_SHUFFLE(1,0,3,2)));
                        sum128 = _mm_add_epi32(sum128, _mm_shuffle_epi32(sum128, _MM_SHUFFLE(2,3,0,1)));
                        total += _mm_extract_epi32(sum128, 0);

                        for (; k < len; k++) {
                            total += temp_b[ny * width + start_x + k];
                        }
                        count += len;
                    }
                    
                    if (count > 0) {
                        b[j * width + i] = (total + count/2) / count;
                    }
                }
            }
        }
    }

    free(temp);
}

int main(int argc, char *argv[]) {
    char *inputs[2] = {NULL, NULL};
    char *output = NULL;
    char *operation = NULL;
    int add_value = 0;
    float contrast_alpha = 1.0;
    int crop_params[4] = {0};
    int rotate_degrees = 0;
    int resize_params[2] = {0};
    float vignette_intensity = 0.0;
    float rgb_scales[3] = {1.0, 1.0, 1.0};
    int blur_params[5] = {0};

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-help") == 0) {
            fprintf(stderr, "\nUsage:\n"
                "  Basic operations:\n"
                "    ./bmpedit -i input.bmp -o output.bmp -op <operation> [parameters]\n\n"
                "  Supported operations:\n"
                "    - add <value>          Adjust brightness (-255 to 255)\n"
                "    - average              Blend two images (requires two -i inputs)\n"
                "    - contrast <alpha>     Adjust contrast (0.1-2.0)\n"
                "    - grayscale            Convert to grayscale\n"
                "    - invert               Invert colors\n"
                "    - crop x y w h         Crop image\n"
                "    - rotate <degrees>     Rotate image (90, -90, 180)\n"
                "    - resize new_w new_h   Resize image\n"
                "    - vignette <intensity> Apply vignette effect (0.0-1.0)\n"
                "    - rgb r_scale g_scale b_scale  Adjust RGB channels\n"
                "    - blur x y w h radius  Apply partial blur\n\n"
                "Examples:\n"
                "  Vignette effect:    ./bmpedit -i input.bmp -o output.bmp -op vignette 0.8\n"
                "  RGB adjustment:     ./bmpedit -i input.bmp -o output.bmp -op rgb 1.2 0.9 1.0\n"
                "  Partial blur:       ./bmpedit -i input.bmp -o output.bmp -op blur 100 100 200 200 3\n");
            return 1;
        }
        else if (strcmp(argv[i], "-i") == 0) {
            if (i+1 >= argc) {
                fprintf(stderr, "Error: Missing argument for -i option\n");
                return 1;
            }
            if (!inputs[0]) {
                inputs[0] = argv[++i];
            }
            else {
                inputs[1] = argv[++i];
            }
        } else if (strcmp(argv[i], "-o") == 0) {
            if (i+1 >= argc) {
                fprintf(stderr, "Error: Missing argument for -o option\n");
                return 1;
            }
            if(output == NULL) {
                output = argv[++i];
            }
            else {
                fprintf(stderr, "Error: Too much output files\n");
                return 1;
            }
        } else if (strcmp(argv[i], "-op") == 0) {
            if (i+1 >= argc) {
                fprintf(stderr, "Error: Missing operation specification\n");
                return 1;
            }
            operation = argv[++i];
            
            if (strcmp(operation, "add") == 0) {
                if (i+1 >= argc) {
                    fprintf(stderr, "Error: Missing value for brightness adjustment\n");
                    return 1;
                }
                add_value = atoi(argv[++i]);
                if (add_value < -255 || add_value > 255) {
                    fprintf(stderr, "Error: Brightness value out of range (-255 to 255)\n");
                    return 1;
                }
            }
            else if (strcmp(operation, "contrast") == 0) {
                if (i+1 >= argc) {
                    fprintf(stderr, "Error: Missing contrast alpha value\n");
                    return 1;
                }
                contrast_alpha = atof(argv[++i]);
                if (contrast_alpha <= 0) {
                    fprintf(stderr, "Error: Contrast alpha must be positive\n");
                    return 1;
                }
            }
            else if (strcmp(operation, "crop") == 0) {
                for (int j = 0; j < 4; j++) {
                    if (i+1 >= argc) {
                        fprintf(stderr, "Error: Insufficient parameters for crop operation\n");
                        return 1;
                    }
                    crop_params[j] = atoi(argv[++i]);
                }
            }
            else if (strcmp(operation, "rotate") == 0) {
                if (i+1 >= argc) {
                    fprintf(stderr, "Error: Missing rotation degrees\n");
                    return 1;
                }
                rotate_degrees = atoi(argv[++i]);
                if (rotate_degrees != 90 && rotate_degrees != -90 && rotate_degrees != 180) {
                    fprintf(stderr, "Error: Invalid rotation degrees. Only 90, -90, and 180 are supported\n");
                    return 1;
                }
            }
            else if (strcmp(operation, "resize") == 0) {
                for (int j = 0; j < 2; j++) {
                    if (i+1 >= argc) {
                        fprintf(stderr, "Error: Insufficient parameters for resize operation\n");
                        return 1;
                    }
                    resize_params[j] = atoi(argv[++i]);
                    if (resize_params[j] <= 0) {
                        fprintf(stderr, "Error: Invalid resize dimensions. Values must be positive\n");
                        return 1;
                    }
                }
            }
            else if (strcmp(operation, "vignette") == 0) {
                if (i+1 >= argc) {
                    fprintf(stderr, "Error: Missing vignette intensity value\n");
                    return 1;
                }
                vignette_intensity = atof(argv[++i]);
                if (vignette_intensity < 0 || vignette_intensity > 1.0) {
                    fprintf(stderr, "Error: Vignette intensity must be between 0.0 and 1.0\n");
                    return 1;
                }
            }
            else if (strcmp(operation, "rgb") == 0) {
                for (int j = 0; j < 3; j++) {
                    if (i+1 >= argc) {
                        fprintf(stderr, "Error: Insufficient parameters for RGB scaling\n");
                        return 1;
                    }
                    rgb_scales[j] = atof(argv[++i]);
                    if (rgb_scales[j] < 0) {
                        fprintf(stderr, "Error: RGB scale values must be non-negative\n");
                        return 1;
                    }
                }
            }
            else if (strcmp(operation, "blur") == 0) {
                for (int j = 0; j < 5; j++) {
                    if (i+1 >= argc) {
                        fprintf(stderr, "Error: Insufficient parameters for blur operation\n");
                        return 1;
                    }
                    blur_params[j] = atoi(argv[++i]);
                }
            }
            else {
                fprintf(stderr, "Error: Unknown operation '%s'\n", operation);
                return 1;
            }
        }
        else {
            fprintf(stderr, "Error: Unknown option '%s'\n", argv[i]);
            return 1;
        }
    }

    if (!inputs[0] || !output || !operation) {
        fprintf(stderr, "\nUsage:\n"
                "  Basic operations:\n"
                "    ./bmpedit -i input.bmp -o output.bmp -op <operation> [parameters]\n\n"
                "  Supported operations:\n"
                "    - add <value>          Adjust brightness (-255 to 255)\n"
                "    - average              Blend two images (requires two -i inputs)\n"
                "    - contrast <alpha>     Adjust contrast (0.1-2.0)\n"
                "    - grayscale            Convert to grayscale\n"
                "    - invert               Invert colors\n"
                "    - crop x y w h         Crop image\n"
                "    - rotate <degrees>     Rotate image (90, -90, 180)\n"
                "    - resize new_w new_h   Resize image\n"
                "    - vignette <intensity> Apply vignette effect (0.0-1.0)\n"
                "    - rgb r_scale g_scale b_scale  Adjust RGB channels\n"
                "    - blur x y w h radius  Apply partial blur\n\n"
                "Examples:\n"
                "  Vignette effect:    ./bmpedit -i input.bmp -o output.bmp -op vignette 0.8\n"
                "  RGB adjustment:     ./bmpedit -i input.bmp -o output.bmp -op rgb 1.2 0.9 1.0\n"
                "  Partial blur:       ./bmpedit -i input.bmp -o output.bmp -op blur 100 100 200 200 3\n");
        return 1;
    }

    char *error = NULL;
    BMPImage *img1 = read_bmp(inputs[0], &error);
    if (!img1) {
        fprintf(stderr, "Error: %s\n", error);
        return 1;
    }

    BMPImage *img2 = NULL;
    BMPImage *result = img1;
    struct timespec start, end;
    if((strcmp(operation, "average") != 0) && inputs[1] != NULL) {
        fprintf(stderr, "Error: Too much inputs\n");
        free_bmp(img1);
        return 1;
    }
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    if (strcmp(operation, "add") == 0) {
        if(add_value < -255 || add_value > 255) {
            fprintf(stderr, "Error: Brightness value out of range (-255 to 255)\n");
            free_bmp(img1);
            return 1;
        }
        adjust_brightness(img1, add_value);
    } else if (strcmp(operation, "average") == 0) {
        if (!inputs[1]) {
            fprintf(stderr, "Error: Blending requires two input files\n");
            free_bmp(img1);
            return 1;
        }
        img2 = read_bmp(inputs[1], &error);
        if (!img2) {
            fprintf(stderr, "Error: %s\n", error);
            free_bmp(img1);
            return 1;
        }
        result = blend_images(img1, img2, &error);
        if (!result) {
            fprintf(stderr, "Error: %s\n", error);
            free_bmp(img1);
            free_bmp(img2);
            return 1;
        }
        free_bmp(img1);
        free_bmp(img2);
    } else if (strcmp(operation, "contrast") == 0) {
        if(contrast_alpha < 0.1 || contrast_alpha > 2.0) {
            fprintf(stderr, "Error: Contrast alpha must be between 0.1 and 2.0\n");
            free_bmp(img1);
            return 1;
        }
        adjust_contrast(img1, contrast_alpha);
    } else if (strcmp(operation, "grayscale") == 0) {
        convert_to_grayscale(img1);
    } else if (strcmp(operation, "invert") == 0) {
        invert_colors(img1);
    } else if (strcmp(operation, "crop") == 0) {
        if(crop_params[0] < 0 || crop_params[1] < 0 ||
           crop_params[2] <= 0 || crop_params[3] <= 0 ||
           crop_params[0] + crop_params[2] > img1->width ||
           crop_params[1] + crop_params[3] > img1->height) {
            fprintf(stderr, "Error: Invalid crop parameters\n");
            free_bmp(img1);
            return 1;
        }
        result = crop_image(img1, crop_params[0], crop_params[1], 
                           crop_params[2], crop_params[3], &error);
        if (!result) {
            fprintf(stderr, "Error: %s\n", error);
            free_bmp(img1);
            return 1;
        }
        free_bmp(img1);
    } else if (strcmp(operation, "rotate") == 0) {
        result = rotate_image(img1, rotate_degrees, &error);
        if (!result) {
            fprintf(stderr, "Error: %s\n", error);
            free_bmp(img1);
            return 1;
        }
        free_bmp(img1);
    } else if (strcmp(operation, "resize") == 0) {
        if(resize_params[0] <= 0 || resize_params[1] <= 0) {
            fprintf(stderr, "Error: Invalid resize dimensions\n");
            free_bmp(img1);
            return 1;
        }
        result = resize_image(img1, resize_params[0], resize_params[1], &error);
        if (!result) {
            fprintf(stderr, "Error: %s\n", error);
            free_bmp(img1);
            return 1;
        }
        free_bmp(img1);
    } else if (strcmp(operation, "vignette") == 0) {
        if(vignette_intensity < 0.0 ) {
            fprintf(stderr, "Error: Vignette intensity must be positive\n");
            free_bmp(img1);
            return 1;
        }
        apply_vignette(img1, vignette_intensity);
    } else if (strcmp(operation, "rgb") == 0) {
        adjust_rgb_channels(img1, rgb_scales[0], rgb_scales[1], rgb_scales[2]);
    } else if (strcmp(operation, "blur") == 0) {
        if(blur_params[0] < 0 || blur_params[1] < 0 ||
           blur_params[2] <= 0 || blur_params[3] <= 0 ||
           blur_params[0] + blur_params[2] > img1->width ||
           blur_params[1] + blur_params[3] > img1->height ||
           blur_params[4] <= 0) {
            fprintf(stderr, "Error: Invalid blur parameters\n");
            free_bmp(img1);
            return 1;
        }
        partial_blur(img1, blur_params[0], blur_params[1], 
                    blur_params[2], blur_params[3], blur_params[4]);
    } else {
        fprintf(stderr, "Error: Unsupported operation '%s'\n", operation);
        free_bmp(img1);
        return 1;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double time_spent = (end.tv_sec - start.tv_sec) * 1000.0 
                      + (end.tv_nsec - start.tv_nsec) / 1000000.0;
    printf("\n====== Performance ======\n");
    printf("Total time: %.3f ms\n", time_spent);

    if (!write_bmp(output, result, &error)) {
        fprintf(stderr, "Error: %s\n", error);
        free_bmp(result);
        return 1;
    }

    free_bmp(result);
    return 0;
}