#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <immintrin.h>
#define ALIGN 32

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

void *aligned_malloc(size_t size, char **error) {
    void *ptr = NULL;
    int ret = posix_memalign(&ptr, ALIGN, size);
    if (ret != 0) {
        *error = "内存对齐分配失败";
        return NULL;
    }
    return ptr;
}

void free_bmp(BMPImage *image) {
    if (image) {
        free(image->data);
        free(image);
    }
}

BMPImage *read_bmp(const char *filename, char **error) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        *error = "无法打开文件";
        return NULL;
    }

    BITMAPFILEHEADER file_header;
    if (fread(&file_header, sizeof(file_header), 1, file) != 1) {
        *error = "读取文件头失败";
        fclose(file);
        return NULL;
    }

    if (file_header.bfType != 0x4D42) {
        *error = "不是BMP文件";
        fclose(file);
        return NULL;
    }

    BITMAPINFOHEADER info_header;
    if (fread(&info_header, sizeof(info_header), 1, file) != 1) {
        *error = "读取信息头失败";
        fclose(file);
        return NULL;
    }

    if (info_header.biBitCount != 24 || info_header.biCompression != 0) {
        *error = "仅支持24位未压缩BMP";
        fclose(file);
        return NULL;
    }

    int height = abs(info_header.biHeight);
    int width = info_header.biWidth;
    size_t row_size = ((width * 3 + 3) / 4) * 4;
    size_t data_size = row_size * height;

    unsigned char *pixel_data = malloc(data_size);
    if (!pixel_data) {
        *error = "内存分配失败";
        fclose(file);
        return NULL;
    }

    fseek(file, file_header.bfOffBits, SEEK_SET);
    if (fread(pixel_data, 1, data_size, file) != data_size) {
        *error = "读取像素数据失败";
        free(pixel_data);
        fclose(file);
        return NULL;
    }
    fclose(file);

    if (info_header.biHeight > 0) {
        unsigned char *temp_row = malloc(row_size);
        for (int y = 0; y < height / 2; y++) {
            unsigned char *row1 = pixel_data + y * row_size;
            unsigned char *row2 = pixel_data + (height - 1 - y) * row_size;
            memcpy(temp_row, row1, row_size);
            memcpy(row1, row2, row_size);
            memcpy(row2, temp_row, row_size);
        }
        free(temp_row);
    }

    int unpadded_row_size = width * 3;
    unsigned char *unpadded_data = malloc(unpadded_row_size * height);
    if (!unpadded_data) {
        *error = "内存分配失败";
        free(pixel_data);
        return NULL;
    }

    for (int y = 0; y < height; y++) {
        memcpy(unpadded_data + y * unpadded_row_size,
               pixel_data + y * row_size,
               unpadded_row_size);
    }
    free(pixel_data);

    BMPImage *image = malloc(sizeof(BMPImage));
    image->file_header = file_header;
    image->info_header = info_header;
    image->width = width;
    image->height = height;
    image->data = unpadded_data;

    return image;
}

int write_bmp(const char *filename, BMPImage *image, char **error) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        *error = "无法创建文件";
        return 0;
    }

    size_t row_size = ((image->width * 3 + 3) / 4) * 4;
    size_t data_size = row_size * image->height;

    unsigned char *padded_data = malloc(data_size);
    if (!padded_data) {
        *error = "内存分配失败";
        fclose(file);
        return 0;
    }

    memset(padded_data, 0, data_size);

    if (image->info_header.biHeight > 0) {
        for (int y = 0; y < image->height; y++) {
            unsigned char *src_row = image->data + (image->height - 1 - y) * image->width * 3;
            memcpy(padded_data + y * row_size, src_row, image->width * 3);
        }
    } else {
        for (int y = 0; y < image->height; y++) {
            memcpy(padded_data + y * row_size,
                   image->data + y * image->width * 3,
                   image->width * 3);
        }
    }

    image->file_header.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + data_size;
    image->info_header.biSizeImage = data_size;
    image->file_header.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
    image->info_header.biHeight = image->height;

    if (fwrite(&image->file_header, sizeof(BITMAPFILEHEADER), 1, file) != 1 ||
        fwrite(&image->info_header, sizeof(BITMAPINFOHEADER), 1, file) != 1 ||
        fwrite(padded_data, 1, data_size, file) != data_size) {
        *error = "写入失败";
        free(padded_data);
        fclose(file);
        return 0;
    }

    free(padded_data);
    fclose(file);
    return 1;
}

// ============== 图像处理函数 ==============
void adjust_brightness(BMPImage *image, int value) {
    unsigned char *p = image->data;
    const int total_bytes = image->width * image->height * 3;
    
    // 限制value范围至[-255, 255]
    value = value > 255 ? 255 : (value < -255 ? -255 : value);

    // 计算对齐的字节数（32字节的整数倍）
    int aligned_bytes = total_bytes - (total_bytes % 32);

    // SIMD+OpenMP优化部分（仅处理对齐块）
    #pragma omp parallel for
    for (int i = 0; i < aligned_bytes; i += 32) {
        // 加载未对齐的256位数据
        __m256i data = _mm256_loadu_si256((__m256i*)(p + i));

        // 将8位数据扩展为16位（零扩展）
        __m256i lo = _mm256_cvtepu8_epi16(_mm256_extracti128_si256(data, 0));
        __m256i hi = _mm256_cvtepu8_epi16(_mm256_extracti128_si256(data, 1));

        // 创建广播的value向量
        __m256i val_vec = _mm256_set1_epi16(value);

        // 执行16位加法
        lo = _mm256_add_epi16(lo, val_vec);
        hi = _mm256_add_epi16(hi, val_vec);

        // 饱和到0-255范围
        lo = _mm256_max_epi16(lo, _mm256_setzero_si256());
        lo = _mm256_min_epi16(lo, _mm256_set1_epi16(255));
        hi = _mm256_max_epi16(hi, _mm256_setzero_si256());
        hi = _mm256_min_epi16(hi, _mm256_set1_epi16(255));

        // 将16位数据重新打包为8位
        __m256i packed = _mm256_packus_epi16(lo, hi);
        packed = _mm256_permute4x64_epi64(packed, _MM_SHUFFLE(3,1,2,0));

        // 存储结果
        _mm256_storeu_si256((__m256i*)(p + i), packed);
    }

    // 处理剩余不足32字节的部分（逐字节处理）
    for (int i = aligned_bytes; i < total_bytes; ++i) {
        int new_val = p[i] + value;
        p[i] = (unsigned char)(new_val > 255 ? 255 : (new_val < 0 ? 0 : new_val));
    }
}

BMPImage *blend_images(BMPImage *img1, BMPImage *img2, char **error) {
    if (img1->width != img2->width || img1->height != img2->height) {
        *error = "图像尺寸不匹配";
        return NULL;
    }

    BMPImage *result = malloc(sizeof(BMPImage));
    if (!result) {
        *error = "内存分配失败";
        return NULL;
    }

    *result = *img1;
    const int total_pixels = img1->width * img1->height;
    const int total_bytes = total_pixels * 3;
    
    result->data = malloc(total_bytes);
    if (!result->data) {
        *error = "内存分配失败";
        free(result);
        return NULL;
    }

    unsigned char *p1 = img1->data;
    unsigned char *p2 = img2->data;
    unsigned char *pr = result->data;

    // SIMD+OpenMP优化部分
    #pragma omp parallel for simd schedule(simd:static) aligned(pr, p1, p2:32)
    for (int i = 0; i < total_bytes; i += 32) {
        // 加载256位数据块
        __m256i data1 = _mm256_loadu_si256((__m256i*)(p1 + i));
        __m256i data2 = _mm256_loadu_si256((__m256i*)(p2 + i));

        // 并行平均计算：使用移位代替除法
        __m256i avg = _mm256_avg_epu8(data1, data2);

        // 存储结果
        _mm256_storeu_si256((__m256i*)(pr + i), avg);
    }

    // 处理剩余不足32字节的部分
    #pragma omp parallel for
    for (int i = (total_bytes / 32) * 32; i < total_bytes; ++i) {
        pr[i] = (p1[i] + p2[i]) / 2;
    }

    return result;
}


void adjust_contrast(BMPImage *image, float alpha) {
    unsigned char *p = image->data;
    const int total_bytes = image->width * image->height * 3;
    
    // 使用Q15定点数格式（32768倍放大）
    const int alpha_q15 = (int)(alpha * 16384);
    const __m256i _128_vec = _mm256_set1_epi16(128);

    // 计算对齐的字节数（32字节的整数倍）
    const int aligned_bytes = total_bytes - (total_bytes % 32);

    // AVX2常数初始化
    const __m256i alpha_vec = _mm256_set1_epi16(alpha_q15);

    #pragma omp parallel for
    for (int i = 0; i < aligned_bytes; i += 32) {
        // 加载未对齐的256位数据
        __m256i data = _mm256_loadu_si256((__m256i*)(p + i));
        // printf("data-0: %d\n", _mm256_extract_epi8(data, 0));
        // 将8位数据扩展为16位（带符号扩展）
        __m256i lo = _mm256_cvtepu8_epi16(_mm256_extracti128_si256(data, 0));
        __m256i hi = _mm256_cvtepu8_epi16(_mm256_extracti128_si256(data, 1));

        // printf("lo-0: %hd\n", _mm256_extract_epi16(lo, 0));
        // printf("hi-0: %hd\n", _mm256_extract_epi16(hi, 0));

        // 减去128中心值（转为有符号运算）
        lo = _mm256_sub_epi16(lo, _128_vec);
        hi = _mm256_sub_epi16(hi, _128_vec);

        // 将数据缩放到Q15范围（左移8位，等效乘256）
        // printf("lo00: %hd\n", _mm256_extract_epi16(lo, 0));
        // printf("hi00: %hd\n", _mm256_extract_epi16(hi, 0));

        lo = _mm256_slli_epi16(lo, 8);  // 左移8位，扩展为Q15
        hi = _mm256_slli_epi16(hi, 8);
        // printf("lo0: %hd\n", _mm256_extract_epi16(lo, 0));
        // printf("hi0: %hd\n", _mm256_extract_epi16(hi, 0));

        // Q15定点数乘法（自动右移15位）
        // printf("lo: %hd\n", _mm256_extract_epi16(lo, 0));
        // printf("hi: %hd\n", _mm256_extract_epi16(hi, 0));
        // printf("alpha: %hd\n", alpha_q15);
        lo = _mm256_mulhrs_epi16(lo, alpha_vec);
        hi = _mm256_mulhrs_epi16(hi, alpha_vec);
        // printf("lo1: %hd\n", _mm256_extract_epi16(lo, 0));
        // printf("hi1: %hd\n", _mm256_extract_epi16(hi, 0));

        // 恢复缩放（右移8位，或通过后续操作隐式处理）
        lo = _mm256_srai_epi16(lo, 8);  // 算术右移8位
        hi = _mm256_srai_epi16(hi, 8);
        // printf("lo2: %hd\n", _mm256_extract_epi16(lo, 0));
        // printf("hi2: %hd\n", _mm256_extract_epi16(hi, 0));

        // 加回128并饱和处理
        lo = _mm256_add_epi16(lo, lo);
        hi = _mm256_add_epi16(hi, hi);
        lo = _mm256_add_epi16(lo, _128_vec);
        hi = _mm256_add_epi16(hi, _128_vec);

        // printf("lo-3: %hd\n", _mm256_extract_epi16(lo, 0));
        // printf("hi-3: %hd\n", _mm256_extract_epi16(hi, 0));

        // 打包为8位无符号数
        __m256i packed = _mm256_packus_epi16(lo, hi);
        packed = _mm256_permute4x64_epi64(packed, _MM_SHUFFLE(3,1,2,0));
        // printf("packed-3: %d\n", _mm256_extract_epi8(packed, 0));
        // 存储结果
        _mm256_storeu_si256((__m256i*)(p + i), packed);
    }

    // 处理剩余不足32字节的部分
    #pragma omp parallel for
    for (int i = aligned_bytes; i < total_bytes; ++i) {
        int val = ((p[i] - 128) * alpha_q15 + 0x4000) >> 15;  // 正确舍入
        val += 128;
        p[i] = (unsigned char)(val > 255 ? 255 : (val < 0 ? 0 : val));
    }
}

void convert_to_grayscale(BMPImage *image) {
    unsigned char *p = image->data;
    for (int i = 0; i < image->width * image->height; i++) {
        unsigned char r = p[2], g = p[1], b = p[0];
        unsigned char gray = 0.299 * r + 0.587 * g + 0.114 * b;
        p[0] = p[1] = p[2] = gray;
        p += 3;
    }
}

void invert_colors(BMPImage *image) {
    unsigned char *p = image->data;
    const int total_bytes = image->width * image->height * 3;
    
    // 计算对齐的字节数（32字节的整数倍）
    const int aligned_bytes = (total_bytes / 32) * 32;

    // SIMD+OpenMP优化部分（处理对齐块）
    #pragma omp parallel for
    for (int i = 0; i < aligned_bytes; i += 32) {
        // 加载未对齐的256位数据
        __m256i data = _mm256_loadu_si256((__m256i*)(p + i));
        
        // 创建全255向量
        const __m256i ones = _mm256_set1_epi8(0xFF);
        
        // 计算反转颜色：255 - data
        const __m256i inverted = _mm256_sub_epi8(ones, data);
        
        // 存储结果
        _mm256_storeu_si256((__m256i*)(p + i), inverted);
    }

    // 处理剩余不足32字节的部分（串行处理）
    for (int i = aligned_bytes; i < total_bytes; ++i) {
        p[i] = 255 - p[i];
    }
}

BMPImage *crop_image(BMPImage *image, int x, int y, int w, int h, char **error) {
    if (x < 0 || y < 0 || x + w > image->width || y + h > image->height) {
        *error = "裁剪区域超出图像范围";
        return NULL;
    }

    BMPImage *result = malloc(sizeof(BMPImage));
    if (!result) {
        *error = "内存分配失败";
        return NULL;
    }

    result->file_header = image->file_header;
    result->info_header = image->info_header;
    result->info_header.biWidth = w;
    result->info_header.biHeight = h;
    result->width = w;
    result->height = h;

    result->data = malloc(w * h * 3);
    if (!result->data) {
        *error = "内存分配失败";
        free(result);
        return NULL;
    }

    for (int row = 0; row < h; row++) {
        unsigned char *src = image->data + ((y + row) * image->width + x) * 3;
        unsigned char *dst = result->data + row * w * 3;
        memcpy(dst, src, w * 3);
    }

    return result;
}

BMPImage *rotate_image(BMPImage *image, int degrees, char **error) {
    if (degrees != 90 && degrees != -90 && degrees != 180) {
        *error = "仅支持90度、-90度和180度旋转";
        return NULL;
    }

    BMPImage *result = malloc(sizeof(BMPImage));
    if (!result) {
        *error = "内存分配失败";
        return NULL;
    }

    result->file_header = image->file_header;
    result->info_header = image->info_header;

    if (degrees == 90 || degrees == -90) {
        result->width = image->height;
        result->height = image->width;
        result->info_header.biWidth = image->height;
        result->info_header.biHeight = image->width;
    } else {
        result->width = image->width;
        result->height = image->height;
    }

    result->data = malloc(result->width * result->height * 3);
    if (!result->data) {
        *error = "内存分配失败";
        free(result);
        return NULL;
    }

    for (int y = 0; y < image->height; y++) {
        for (int x = 0; x < image->width; x++) {
            unsigned char *src = image->data + (y * image->width + x) * 3;
            int dst_x, dst_y;
            if (degrees == 90) {
                dst_x = image->height - 1 - y;
                dst_y = x;
            } else if (degrees == -90) {
                dst_x = y;
                dst_y = image->width - 1 - x;
            } else {
                dst_x = image->width - 1 - x;
                dst_y = image->height - 1 - y;
            }
            unsigned char *dst = result->data + (dst_y * result->width + dst_x) * 3;
            memcpy(dst, src, 3);
        }
    }

    return result;
}

BMPImage *resize_image(BMPImage *image, int new_width, int new_height, char **error) {
    if (new_width <= 0 || new_height <= 0) {
        *error = "缩放尺寸必须大于零";
        return NULL;
    }

    BMPImage *result = malloc(sizeof(BMPImage));
    if (!result) {
        *error = "内存分配失败";
        return NULL;
    }

    result->file_header = image->file_header;
    result->info_header = image->info_header;
    result->width = new_width;
    result->height = new_height;
    result->info_header.biWidth = new_width;
    result->info_header.biHeight = new_height;

    result->data = malloc(new_width * new_height * 3);
    if (!result->data) {
        *error = "内存分配失败";
        free(result);
        return NULL;
    }

    float x_ratio = (float)image->width / new_width;
    float y_ratio = (float)image->height / new_height;

    for (int y = 0; y < new_height; y++) {
        for (int x = 0; x < new_width; x++) {
            int src_x = (int)(x * x_ratio);
            int src_y = (int)(y * y_ratio);
            src_x = (src_x >= image->width) ? image->width - 1 : src_x;
            src_y = (src_y >= image->height) ? image->height - 1 : src_y;

            unsigned char *src_pixel = image->data + (src_y * image->width + src_x) * 3;
            unsigned char *dst_pixel = result->data + (y * new_width + x) * 3;
            memcpy(dst_pixel, src_pixel, 3);
        }
    }

    return result;
}

// ============== 新增功能 ==============
void apply_vignette(BMPImage *image, float intensity) {
    int center_x = image->width / 2;
    int center_y = image->height / 2;
    float max_dist = sqrtf(center_x * center_x + center_y * center_y);

    unsigned char *p = image->data;
    for (int y = 0; y < image->height; y++) {
        for (int x = 0; x < image->width; x++) {
            float dx = x - center_x;
            float dy = y - center_y;
            float distance = sqrtf(dx*dx + dy*dy) / max_dist;
            float factor = 1.0 - intensity * distance * distance;

            for (int c = 0; c < 3; c++) {
                int new_val = p[c] * factor;
                p[c] = (unsigned char)fmax(fmin(new_val, 255), 0);
            }
            p += 3;
        }
    }
}

void adjust_rgb_channels(BMPImage *image, float r_scale, float g_scale, float b_scale) {
    unsigned char *p = image->data;
    for (int i = 0; i < image->width * image->height * 3; i += 3) {
        p[i+2] = fmin(p[i+2] * r_scale, 255); // Red
        p[i+1] = fmin(p[i+1] * g_scale, 255); // Green
        p[i]   = fmin(p[i]   * b_scale, 255); // Blue
    }
}

void partial_blur(BMPImage *image, int x, int y, int w, int h, int radius) {
    unsigned char *temp = malloc(image->width * image->height * 3);
    memcpy(temp, image->data, image->width * image->height * 3);

    for (int j = y; j < y + h; j++) {
        for (int i = x; i < x + w; i++) {
            int sum_r = 0, sum_g = 0, sum_b = 0, count = 0;
            
            for (int dy = -radius; dy <= radius; dy++) {
                for (int dx = -radius; dx <= radius; dx++) {
                    int nx = i + dx;
                    int ny = j + dy;
                    if (nx >= 0 && nx < image->width && ny >= 0 && ny < image->height) {
                        unsigned char *p = temp + (ny * image->width + nx) * 3;
                        sum_r += p[2];
                        sum_g += p[1];
                        sum_b += p[0];
                        count++;
                    }
                }
            }
            
            if (count > 0) {
                unsigned char *dst = image->data + (j * image->width + i) * 3;
                dst[2] = sum_r / count;
                dst[1] = sum_g / count;
                dst[0] = sum_b / count;
            }
        }
    }
    free(temp);
}

// ============== 主函数 ==============
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
        if (strcmp(argv[i], "-i") == 0) {
            if (!inputs[0]) inputs[0] = argv[++i];
            else inputs[1] = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0) {
            output = argv[++i];
        } else if (strcmp(argv[i], "-op") == 0) {
            operation = argv[++i];
            if (strcmp(operation, "add") == 0 && i+1 < argc) {
                add_value = atoi(argv[++i]);
            } else if (strcmp(operation, "contrast") == 0 && i+1 < argc) {
                contrast_alpha = atof(argv[++i]);
            } else if (strcmp(operation, "crop") == 0) {
                for (int j = 0; j < 4 && i+1 < argc; j++) {
                    crop_params[j] = atoi(argv[++i]);
                }
            } else if (strcmp(operation, "rotate") == 0 && i+1 < argc) {
                rotate_degrees = atoi(argv[++i]);
            } else if (strcmp(operation, "resize") == 0) {
                for (int j = 0; j < 2 && i+1 < argc; j++) {
                    resize_params[j] = atoi(argv[++i]);
                }
            } else if (strcmp(operation, "vignette") == 0 && i+1 < argc) {
                vignette_intensity = atof(argv[++i]);
            } else if (strcmp(operation, "rgb") == 0) {
                for (int j = 0; j < 3 && i+1 < argc; j++) {
                    rgb_scales[j] = atof(argv[++i]);
                }
            } else if (strcmp(operation, "blur") == 0) {
                for (int j = 0; j < 5 && i+1 < argc; j++) {
                    blur_params[j] = atoi(argv[++i]);
                }
            }
        }
    }

    if (!inputs[0] || !output || !operation) {
        fprintf(stderr, "参数错误\n用法示例：\n"
                        "  - 渐晕效果:    ./bmpedit -i input.bmp -o output.bmp -op vignette 0.8\n"
                        "  - RGB通道调整: ./bmpedit -i input.bmp -o output.bmp -op rgb 1.2 0.9 1.0\n"
                        "  - 区域模糊:    ./bmpedit -i input.bmp -o output.bmp -op blur 100 100 200 200 3\n"
                        "  - 其他原有操作参数保持不变...\n");
        return 1;
    }

    char *error = NULL;
    BMPImage *img1 = read_bmp(inputs[0], &error);
    if (!img1) {
        fprintf(stderr, "错误: %s\n", error);
        return 1;
    }

    BMPImage *img2 = NULL;
    BMPImage *result = img1;

    if (strcmp(operation, "add") == 0) {
        adjust_brightness(img1, add_value);
    } else if (strcmp(operation, "average") == 0) {
        if (!inputs[1]) {
            fprintf(stderr, "混合需要两个输入文件\n");
            free_bmp(img1);
            return 1;
        }
        img2 = read_bmp(inputs[1], &error);
        if (!img2) {
            fprintf(stderr, "错误: %s\n", error);
            free_bmp(img1);
            return 1;
        }
        result = blend_images(img1, img2, &error);
        if (!result) {
            fprintf(stderr, "错误: %s\n", error);
            free_bmp(img1);
            free_bmp(img2);
            return 1;
        }
        free_bmp(img1);
        free_bmp(img2);
    } else if (strcmp(operation, "contrast") == 0) {
        adjust_contrast(img1, contrast_alpha);
    } else if (strcmp(operation, "grayscale") == 0) {
        convert_to_grayscale(img1);
    } else if (strcmp(operation, "invert") == 0) {
        invert_colors(img1);
    } else if (strcmp(operation, "crop") == 0) {
        result = crop_image(img1, crop_params[0], crop_params[1], 
                           crop_params[2], crop_params[3], &error);
        if (!result) {
            fprintf(stderr, "错误: %s\n", error);
            free_bmp(img1);
            return 1;
        }
        free_bmp(img1);
    } else if (strcmp(operation, "rotate") == 0) {
        result = rotate_image(img1, rotate_degrees, &error);
        if (!result) {
            fprintf(stderr, "错误: %s\n", error);
            free_bmp(img1);
            return 1;
        }
        free_bmp(img1);
    } else if (strcmp(operation, "resize") == 0) {
        if (resize_params[0] <= 0 || resize_params[1] <= 0) {
            fprintf(stderr, "错误: 缩放尺寸必须大于零\n");
            free_bmp(img1);
            return 1;
        }
        result = resize_image(img1, resize_params[0], resize_params[1], &error);
        if (!result) {
            fprintf(stderr, "错误: %s\n", error);
            free_bmp(img1);
            return 1;
        }
        free_bmp(img1);
    } else if (strcmp(operation, "vignette") == 0) {
        apply_vignette(img1, vignette_intensity);
    } else if (strcmp(operation, "rgb") == 0) {
        adjust_rgb_channels(img1, rgb_scales[0], rgb_scales[1], rgb_scales[2]);
    } else if (strcmp(operation, "blur") == 0) {
        partial_blur(img1, blur_params[0], blur_params[1], 
                    blur_params[2], blur_params[3], blur_params[4]);
    } else {
        fprintf(stderr, "不支持的操作: %s\n", operation);
        free_bmp(img1);
        return 1;
    }

    if (!write_bmp(output, result, &error)) {
        fprintf(stderr, "错误: %s\n", error);
        free_bmp(result);
        return 1;
    }

    free_bmp(result);
    return 0;
}