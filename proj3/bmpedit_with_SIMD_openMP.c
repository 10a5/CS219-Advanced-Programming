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
    unsigned char *data;  // 平面存储: R平面 | G平面 | B平面
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
        *error = "无法打开文件";
        return NULL;
    }

    // 读取文件头
    BITMAPFILEHEADER file_header;
    if (fread(&file_header, sizeof(file_header), 1, file) != 1) {
        *error = "读取文件头失败";
        fclose(file);
        return NULL;
    }

    // 验证BMP格式
    if (file_header.bfType != 0x4D42) {
        *error = "不是BMP文件";
        fclose(file);
        return NULL;
    }

    // 读取信息头
    BITMAPINFOHEADER info_header;
    if (fread(&info_header, sizeof(info_header), 1, file) != 1) {
        *error = "读取信息头失败";
        fclose(file);
        return NULL;
    }

    // 验证格式支持
    if (info_header.biBitCount != 24 || info_header.biCompression != 0) {
        *error = "仅支持24位未压缩BMP";
        fclose(file);
        return NULL;
    }

    // 计算实际尺寸
    const int src_height = abs(info_header.biHeight);
    const int width = info_header.biWidth;
    const size_t row_size = ((width * 3 + 3) / 4) * 4;
    const size_t data_size = row_size * src_height;

    // 读取原始像素数据
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

    // 创建平面存储空间
    const size_t plane_size = width * src_height;
    unsigned char *planar_data = malloc(plane_size * 3);
    if (!planar_data) {
        *error = "内存分配失败";
        free(pixel_data);
        return NULL;
    }

    // 设置平面指针
    unsigned char *r_plane = planar_data;
    unsigned char *g_plane = planar_data + plane_size;
    unsigned char *b_plane = planar_data + plane_size * 2;

    // 转换到平面格式并处理方向
    if (info_header.biHeight > 0) {
        // 自下而上的原始数据需要翻转
        for (int y = 0; y < src_height; y++) {
            const unsigned char *src_row = pixel_data + (src_height - 1 - y) * row_size;
            for (int x = 0; x < width; x++) {
                const int dst_idx = y * width + x;
                b_plane[dst_idx] = src_row[x * 3];     // B
                g_plane[dst_idx] = src_row[x * 3 + 1]; // G
                r_plane[dst_idx] = src_row[x * 3 + 2]; // R
            }
        }
    } else {
        // 自上而下的原始数据直接复制
        for (int y = 0; y < src_height; y++) {
            const unsigned char *src_row = pixel_data + y * row_size;
            for (int x = 0; x < width; x++) {
                const int dst_idx = y * width + x;
                b_plane[dst_idx] = src_row[x * 3];     // B
                g_plane[dst_idx] = src_row[x * 3 + 1]; // G
                r_plane[dst_idx] = src_row[x * 3 + 2]; // R
            }
        }
    }
    free(pixel_data);

    // 创建图像对象
    BMPImage *image = malloc(sizeof(BMPImage));
    image->file_header = file_header;
    image->info_header = info_header;
    image->width = width;
    image->height = src_height;  // 存储为正值
    image->data = planar_data;

    // 强制头信息高度为负值（表示自上而下存储）
    image->info_header.biHeight = -src_height;
    return image;
}

int write_bmp(const char *filename, BMPImage *image, char **error) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        *error = "无法创建文件";
        return 0;
    }

    // 准备交错格式数据
    const size_t row_size = ((image->width * 3 + 3) / 4) * 4;
    const size_t data_size = row_size * image->height;
    unsigned char *padded_data = calloc(1, data_size);
    if (!padded_data) {
        *error = "内存分配失败";
        fclose(file);
        return 0;
    }

    // 获取平面数据指针
    const size_t plane_size = image->width * image->height;
    const unsigned char *r_plane = image->data;
    const unsigned char *g_plane = image->data + plane_size;
    const unsigned char *b_plane = image->data + plane_size * 2;

    // 转换为BGR交错格式（始终自上而下）
    for (int y = 0; y < image->height; y++) {
        unsigned char *dst_row = padded_data + y * row_size;
        for (int x = 0; x < image->width; x++) {
            const int src_idx = y * image->width + x;
            dst_row[x * 3] = b_plane[src_idx];     // B
            dst_row[x * 3 + 1] = g_plane[src_idx]; // G
            dst_row[x * 3 + 2] = r_plane[src_idx]; // R
        }
    }

    // 更新头信息
    image->file_header.bfSize = sizeof(BITMAPFILEHEADER) + 
                               sizeof(BITMAPINFOHEADER) + 
                               data_size;
    image->info_header.biSizeImage = data_size;
    image->file_header.bfOffBits = sizeof(BITMAPFILEHEADER) + 
                                  sizeof(BITMAPINFOHEADER);
    
    // 保持负高度值表示自上而下存储
    image->info_header.biHeight = -image->height;

    // 写入文件
    int success = 1;
    if (fwrite(&image->file_header, sizeof(BITMAPFILEHEADER), 1, file) != 1 ||
        fwrite(&image->info_header, sizeof(BITMAPINFOHEADER), 1, file) != 1 ||
        fwrite(padded_data, 1, data_size, file) != data_size) {
        *error = "写入失败";
        success = 0;
    }

    free(padded_data);
    fclose(file);
    return success;
}

// ============== 图像处理函数 (修改为平面处理) ==============
void adjust_brightness(BMPImage *image, int value) {
    const size_t total = 3 * image->width * image->height;
    unsigned char *data = image->data;
    const int vector_size = 16;  // 移出并行区域
    const size_t end = total - vector_size;

    #pragma omp parallel
    {
        #pragma omp for
        for (size_t i = 0; i <= end; i += vector_size) {
            // 加载16字节
            __m128i vec = _mm_loadu_si128((__m128i*)(data + i));
            
            // 扩展为两个16位向量
            __m128i low = _mm_cvtepu8_epi16(vec);
            __m128i high = _mm_cvtepu8_epi16(_mm_srli_si128(vec, 8));
            
            // 创建值的向量
            __m128i val_vec = _mm_set1_epi16(value);
            
            // 加法运算
            low = _mm_add_epi16(low, val_vec);
            high = _mm_add_epi16(high, val_vec);
            
            // 饱和运算
            low = _mm_min_epi16(_mm_max_epi16(low, _mm_setzero_si128()), _mm_set1_epi16(255));
            high = _mm_min_epi16(_mm_max_epi16(high, _mm_setzero_si128()), _mm_set1_epi16(255));
            
            // 重新打包为8位
            __m128i res = _mm_packus_epi16(low, high);
            
            // 存储结果
            _mm_storeu_si128((__m128i*)(data + i), res);
        }

        // 处理剩余字节（标量处理）
        #pragma omp for
        for (size_t i = ((total / vector_size) * vector_size); i < total; ++i) {
            int temp = data[i] + value;
            data[i] = (unsigned char)(temp > 255 ? 255 : (temp < 0 ? 0 : temp));
        }
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
    const size_t plane_size = image->width * image->height;
    const size_t total_pixels = plane_size;
    unsigned char *r = image->data;
    unsigned char *g = r + plane_size;
    unsigned char *b = g + plane_size;

    // 灰度系数（IEEE标准）
    const __m256 r_coeff = _mm256_set1_ps(0.299f);
    const __m256 g_coeff = _mm256_set1_ps(0.587f);
    const __m256 b_coeff = _mm256_set1_ps(0.114f);
    const __m256 max_val = _mm256_set1_ps(255.0f);
    const __m256 zero = _mm256_setzero_ps();

    // 并行处理所有像素
    #pragma omp parallel for schedule(guided)
    for (size_t i = 0; i < total_pixels; i += 8) {
        // 加载8个像素的三通道数据
        __m128i r8 = _mm_loadl_epi64((__m128i*)(r + i));
        __m128i g8 = _mm_loadl_epi64((__m128i*)(g + i));
        __m128i b8 = _mm_loadl_epi64((__m128i*)(b + i));

        // 转换为32位整数
        __m256i r32 = _mm256_cvtepu8_epi32(r8);
        __m256i g32 = _mm256_cvtepu8_epi32(g8);
        __m256i b32 = _mm256_cvtepu8_epi32(b8);

        // 转换为浮点数
        __m256 rf = _mm256_cvtepi32_ps(r32);
        __m256 gf = _mm256_cvtepi32_ps(g32);
        __m256 bf = _mm256_cvtepi32_ps(b32);

        // 计算灰度值：gray = R*0.299 + G*0.587 + B*0.114
        __m256 gray_f = _mm256_add_ps(
            _mm256_mul_ps(rf, r_coeff),
            _mm256_mul_ps(gf, g_coeff));
        gray_f = _mm256_add_ps(gray_f, _mm256_mul_ps(bf, b_coeff));

        // 限制到[0, 255]范围
        gray_f = _mm256_min_ps(_mm256_max_ps(gray_f, zero), max_val);

        // 转换回整数
        __m256i gray32 = _mm256_cvtps_epi32(gray_f);
        __m128i gray16 = _mm_packus_epi32(
            _mm256_extractf128_si256(gray32, 0),
            _mm256_extractf128_si256(gray32, 1));
        __m128i gray8 = _mm_packus_epi16(gray16, gray16);

        // 存储到三个通道
        _mm_storel_epi64((__m128i*)(r + i), gray8);
        _mm_storel_epi64((__m128i*)(g + i), gray8);
        _mm_storel_epi64((__m128i*)(b + i), gray8);
    }

    // 处理剩余像素（标量处理）
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
    const int vector_size = 16;  // 移出并行区域
    const size_t end = total - vector_size;

    #pragma omp parallel
    {
        #pragma omp for
        for (size_t i = 0; i <= end; i += vector_size) {
            // 加载16字节
            __m128i vec = _mm_loadu_si128((__m128i*)(data + i));
            
            // 异或运算（按位取反）
            __m128i res = _mm_xor_si128(vec, ones);
            
            // 存储结果
            _mm_storeu_si128((__m128i*)(data + i), res);
        }

        // 处理剩余字节（标量处理）
        #pragma omp for
        for (size_t i = ((total / vector_size) * vector_size); i < total; ++i) {
            data[i] = 255 - data[i];
        }
    }
}

BMPImage *crop_image(BMPImage *image, int x, int y, int w, int h, char **error) {
    if (x < 0 || y < 0 || x + w > image->width || y + h > image->height) {
        *error = "裁剪区域超出图像范围";
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
        *error = "内存分配失败";
        free(result);
        return NULL;
    }

    unsigned char *src_r = image->data;
    unsigned char *src_g = src_r + image->width * image->height;
    unsigned char *src_b = src_g + image->width * image->height;

    unsigned char *dst_r = result->data;
    unsigned char *dst_g = dst_r + new_plane_size;
    unsigned char *dst_b = dst_g + new_plane_size;

    // 使用并行化行处理
    #pragma omp parallel for schedule(guided)
    for (int row = 0; row < h; row++) {
        const int src_y = y + row;
        if (src_y >= 0 && src_y < image->height) {
            for (int col = 0; col < w; col++) {
                const int src_x = x + col;
                if (src_x >= 0 && src_x < image->width) {
                    const int src_idx = src_y * image->width + src_x;
                    const int dst_idx = row * w + col;
                    
                    // 批量内存拷贝优化
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
        *error = "仅支持90度、-90度和180度旋转";
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
        *error = "内存分配失败";
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

    // 分情况处理不同旋转角度
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
// ==================== 优化后的resize函数 ====================
BMPImage *resize_image(BMPImage *image, int new_width, int new_height, char **error) {
    if (new_width <= 0 || new_height <= 0) {
        *error = "缩放尺寸必须大于零";
        return NULL;
    }

    BMPImage *result = malloc(sizeof(BMPImage));
    *result = *image;
    result->width = new_width;
    result->height = new_height;
    result->info_header.biWidth = new_width;
    result->info_header.biHeight = new_height;

    const size_t new_plane_size = new_width * new_height;
    unsigned char *planar_data = aligned_alloc(64, new_plane_size * 3); // 64字节对齐
    if (!planar_data) {
        *error = "内存分配失败";
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

    // 并行化Y轴循环
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

// ============== 新增功能 ==============
void apply_vignette(BMPImage *image, float intensity) {
    const int width = image->width;
    const int height = image->height;
    const size_t plane_size = width * height;
    unsigned char *r = image->data;
    unsigned char *g = image->data + plane_size;
    unsigned char *b = image->data + plane_size * 2;

    // 预计算常量
    const float center_x = (width - 1) * 0.5f;
    const float center_y = (height - 1) * 0.5f;
    const float max_radius_sq = center_x*center_x + center_y*center_y;
    const __m256 intensity_vec = _mm256_set1_ps(intensity / max_radius_sq);
    const __m256 one_vec = _mm256_set1_ps(1.0f);
    const __m256 zero_vec = _mm256_setzero_ps();
    const __m256i mask = _mm256_setr_epi32(0,1,2,3,4,5,6,7);

    // OpenMP并行处理行
    #pragma omp parallel for schedule(guided)
    for (int y = 0; y < height; ++y) {
        const __m256 dy = _mm256_set1_ps(y - center_y);
        const __m256 dy_sq = _mm256_mul_ps(dy, dy);
        unsigned char *r_row = r + y*width;
        unsigned char *g_row = g + y*width;
        unsigned char *b_row = b + y*width;

        // SIMD处理8像素块
        int x = 0;
        for (; x <= width-8; x += 8) {
            // 生成x坐标：x+0到x+7
            __m256i x_idx = _mm256_add_epi32(_mm256_set1_epi32(x), mask);
            __m256 x_f = _mm256_cvtepi32_ps(x_idx);
            
            // 计算dx = x_f - center_x
            __m256 dx = _mm256_sub_ps(x_f, _mm256_set1_ps(center_x));
            __m256 dx_sq = _mm256_mul_ps(dx, dx);
            
            // distance_sq = (dx² + dy²)
            __m256 distance_sq = _mm256_add_ps(dx_sq, dy_sq);
            
            // factor = 1 - intensity * distance_sq / max_radius_sq
            __m256 factor = _mm256_mul_ps(distance_sq, intensity_vec);
            factor = _mm256_sub_ps(one_vec, factor);

            // 加载并处理RGB数据
            __m256i r_data = _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i*)(r_row+x)));
            __m256i g_data = _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i*)(g_row+x)));
            __m256i b_data = _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i*)(b_row+x)));

            __m256 r_f = _mm256_cvtepi32_ps(r_data);
            __m256 g_f = _mm256_cvtepi32_ps(g_data);
            __m256 b_f = _mm256_cvtepi32_ps(b_data);

            // 应用渐晕因子
            r_f = _mm256_mul_ps(r_f, factor);
            g_f = _mm256_mul_ps(g_f, factor);
            b_f = _mm256_mul_ps(b_f, factor);

            // 饱和处理
            r_f = _mm256_max_ps(_mm256_min_ps(r_f, _mm256_set1_ps(255.0f)), zero_vec);
            g_f = _mm256_max_ps(_mm256_min_ps(g_f, _mm256_set1_ps(255.0f)), zero_vec);
            b_f = _mm256_max_ps(_mm256_min_ps(b_f, _mm256_set1_ps(255.0f)), zero_vec);

            // 转换回整数并打包存储
            __m256i r_result = _mm256_cvtps_epi32(r_f);
            __m256i g_result = _mm256_cvtps_epi32(g_f);
            __m256i b_result = _mm256_cvtps_epi32(b_f);
            // printf("r_result-0: %d\n", _mm256_extract_epi32(r_result, 0));
            // printf("g_result-0: %d\n", _mm256_extract_epi32(g_result, 0));
            // printf("b_result-0: %d\n", _mm256_extract_epi32(b_result, 0));
            // printf("r_result-1: %d\n", _mm256_extract_epi32(r_result, 1));
            // printf("g_result-1: %d\n", _mm256_extract_epi32(g_result, 1));
            // printf("b_result-1: %d\n", _mm256_extract_epi32(b_result, 1));
            // printf("r_result-2: %d\n", _mm256_extract_epi32(r_result, 2));
            // printf("g_result-2: %d\n", _mm256_extract_epi32(g_result, 2));
            // printf("b_result-2: %d\n", _mm256_extract_epi32(b_result, 2));
            // printf("r_result-3: %d\n", _mm256_extract_epi32(r_result, 3));
            // printf("g_result-3: %d\n", _mm256_extract_epi32(g_result, 3));
            // printf("b_result-3: %d\n", _mm256_extract_epi32(b_result, 3));

            __m128i low_r = _mm256_castsi256_si128(r_result);
            __m128i high_r = _mm256_extracti128_si256(r_result, 1);
            __m128i low_g = _mm256_castsi256_si128(g_result);
            __m128i high_g = _mm256_extracti128_si256(g_result, 1);
            __m128i low_b = _mm256_castsi256_si128(b_result);
            __m128i high_b = _mm256_extracti128_si256(b_result, 1);

            // 定义洗牌掩码：提取非零的 16 位元素（跳过每个 0）
            const __m128i shuffle_mask = _mm_setr_epi8(
                0, 1, 4, 5, 8, 9, 12, 13, // 提取低 128 位的 a, b, c, d（8 字节）
                0, 1, 4, 5, 8, 9, 12, 13  // 提取高 128 位的 e, f, g, h（8 字节）
            );

            // 对低 128 位和高 128 位进行洗牌操作
            __m128i packed_low_r = _mm_shuffle_epi8(low_r, shuffle_mask);  // 低 64 位：a, b, c, d
            __m128i packed_high_r = _mm_shuffle_epi8(high_r, shuffle_mask);// 低 64 位：e, f, g, h
            __m128i packed_low_g = _mm_shuffle_epi8(low_g, shuffle_mask);  // 低 64 位：a, b, c, d
            __m128i packed_high_g = _mm_shuffle_epi8(high_g, shuffle_mask);// 低 64 位：e, f, g, h
            __m128i packed_low_b = _mm_shuffle_epi8(low_b, shuffle_mask);  // 低 64 位：a, b, c, d
            __m128i packed_high_b = _mm_shuffle_epi8(high_b, shuffle_mask);// 低 64 位：e, f, g, h

            // printf("packed_low_r-0: %d\n", _mm_extract_epi16(packed_low_r, 0));
            // printf("packed_high_r-0: %d\n", _mm_extract_epi16(packed_high_r, 0));
            // printf("packed_low_g-0: %d\n", _mm_extract_epi16(packed_low_g, 0));
            // printf("packed_high_g-0: %d\n", _mm_extract_epi16(packed_high_g, 0));
            // printf("packed_low_b-0: %d\n", _mm_extract_epi16(packed_low_b, 0));
            // printf("packed_high_b-0: %d\n", _mm_extract_epi16(packed_high_b, 0));
            // printf("packed_low_r-1: %d\n", _mm_extract_epi16(packed_low_r, 1));
            // printf("packed_high_r-1: %d\n", _mm_extract_epi16(packed_high_r, 1));
            // printf("packed_low_g-1: %d\n", _mm_extract_epi16(packed_low_g, 1));
            // printf("packed_high_g-1: %d\n", _mm_extract_epi16(packed_high_g, 1));
            // printf("packed_low_b-1: %d\n", _mm_extract_epi16(packed_low_b, 1));
            // printf("packed_high_b-1: %d\n", _mm_extract_epi16(packed_high_b, 1));
            // printf("packed_low_r-2: %d\n", _mm_extract_epi16(packed_low_r, 2));
            // printf("packed_high_r-2: %d\n", _mm_extract_epi16(packed_high_r, 2));
            // printf("packed_low_g-2: %d\n", _mm_extract_epi16(packed_low_g, 2));
            // printf("packed_high_g-2: %d\n", _mm_extract_epi16(packed_high_g, 2));
            // printf("packed_low_b-2: %d\n", _mm_extract_epi16(packed_low_b, 2));
            // printf("packed_high_b-2: %d\n", _mm_extract_epi16(packed_high_b, 2));
            // printf("packed_low_r-3: %d\n", _mm_extract_epi16(packed_low_r, 3));
            // printf("packed_high_r-3: %d\n", _mm_extract_epi16(packed_high_r, 3));
            // printf("packed_low_g-3: %d\n", _mm_extract_epi16(packed_low_g, 3));
            // printf("packed_high_g-3: %d\n", _mm_extract_epi16(packed_high_g, 3));
            // printf("packed_low_b-3: %d\n", _mm_extract_epi16(packed_low_b, 3));
            // printf("packed_high_b-3: %d\n", _mm_extract_epi16(packed_high_b, 3));

            // 合并低 64 位和高 64 位到 128 位中
            __m128i result_r = _mm_unpacklo_epi64(packed_low_r, packed_high_r);
            __m128i result_g = _mm_unpacklo_epi64(packed_low_g, packed_high_g);
            __m128i result_b = _mm_unpacklo_epi64(packed_low_b, packed_high_b);

            // printf("result_r-0: %d\n", _mm_extract_epi16(result_r, 0));
            // printf("result_g-0: %d\n", _mm_extract_epi16(result_g, 0));
            // printf("result_b-0: %d\n", _mm_extract_epi16(result_b, 0));
            // printf("result_r-1: %d\n", _mm_extract_epi16(result_r, 1));
            // printf("result_g-1: %d\n", _mm_extract_epi16(result_g, 1));
            // printf("result_b-1: %d\n", _mm_extract_epi16(result_b, 1));
            // printf("result_r-2: %d\n", _mm_extract_epi16(result_r, 2));
            // printf("result_g-2: %d\n", _mm_extract_epi16(result_g, 2));
            // printf("result_b-2: %d\n", _mm_extract_epi16(result_b, 2));
            // printf("result_r-3: %d\n", _mm_extract_epi16(result_r, 3));
            // printf("result_g-3: %d\n", _mm_extract_epi16(result_g, 3));
            // printf("result_b-3: %d\n", _mm_extract_epi16(result_b, 3));
            // printf("result_r-4: %d\n", _mm_extract_epi16(result_r, 4));
            // printf("result_g-4: %d\n", _mm_extract_epi16(result_g, 4));
            // printf("result_b-4: %d\n", _mm_extract_epi16(result_b, 4));
            // printf("result_r-5: %d\n", _mm_extract_epi16(result_r, 5));
            // printf("result_g-5: %d\n", _mm_extract_epi16(result_g, 5));
            // printf("result_b-5: %d\n", _mm_extract_epi16(result_b, 5));
            // printf("result_r-6: %d\n", _mm_extract_epi16(result_r, 6));
            // printf("result_g-6: %d\n", _mm_extract_epi16(result_g, 6));
            // printf("result_b-6: %d\n", _mm_extract_epi16(result_b, 6));
            // printf("result_r-7: %d\n", _mm_extract_epi16(result_r, 7));
            // printf("result_g-7: %d\n", _mm_extract_epi16(result_g, 7));
            // printf("result_b-7: %d\n", _mm_extract_epi16(result_b, 7));

            const __m128i shuffle_mask1 = _mm_setr_epi8(
                0, 2, 4, 6, 8, 10, 12, 14, // 提取 a(0), b(2), c(4), d(6), e(8), f(10), g(12), h(14)
                -1, -1, -1, -1, -1, -1, -1, -1  // 高8字节填充0（0xFF对应清零）
            );
            __m128i r_r = _mm_shuffle_epi8(result_r, shuffle_mask1);
            __m128i g_r = _mm_shuffle_epi8(result_g, shuffle_mask1);
            __m128i b_r = _mm_shuffle_epi8(result_b, shuffle_mask1);

            _mm_storel_epi64((__m128i*)(r_row+x), r_r);
            _mm_storel_epi64((__m128i*)(g_row+x), g_r);
            _mm_storel_epi64((__m128i*)(b_row+x), b_r);
        }

        // 处理剩余像素（标量）
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

// ==================== RGB通道调整优化（修复版） ====================
void adjust_rgb_channels(BMPImage *image, float r_scale, float g_scale, float b_scale) {
    const size_t plane_size = image->width * image->height;
    const __m256 r_scale_vec = _mm256_set1_ps(r_scale);
    const __m256 g_scale_vec = _mm256_set1_ps(g_scale);
    const __m256 b_scale_vec = _mm256_set1_ps(b_scale);
    const __m256 max_val = _mm256_set1_ps(255.0f);
    const __m256 min_val = _mm256_setzero_ps();

    #pragma omp parallel sections
    {
        // R通道处理
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

        // G通道处理（完整实现）
        #pragma omp section
        {
            unsigned char *g = image->data + plane_size;
            #pragma omp parallel for simd
            for (size_t i = 0; i < plane_size; i += 8) {
                __m128i vec8 = _mm_loadl_epi64((__m128i*)(g + i));
                __m256i vec32 = _mm256_cvtepu8_epi32(vec8);
                __m256 fvec = _mm256_cvtepi32_ps(vec32);
                fvec = _mm256_mul_ps(fvec, g_scale_vec);  // 使用g_scale_vec
                fvec = _mm256_min_ps(_mm256_max_ps(fvec, min_val), max_val);
                __m256i res32 = _mm256_cvtps_epi32(fvec);
                __m128i res16 = _mm_packus_epi32(_mm256_extracti128_si256(res32, 0),
                                               _mm256_extracti128_si256(res32, 1));
                __m128i res8 = _mm_packus_epi16(res16, res16);
                _mm_storel_epi64((__m128i*)(g + i), res8);
            }
        }

        // B通道处理（完整实现）
        #pragma omp section
        {
            unsigned char *b = image->data + plane_size * 2;
            #pragma omp parallel for simd
            for (size_t i = 0; i < plane_size; i += 8) {
                __m128i vec8 = _mm_loadl_epi64((__m128i*)(b + i));
                __m256i vec32 = _mm256_cvtepu8_epi32(vec8);
                __m256 fvec = _mm256_cvtepi32_ps(vec32);
                fvec = _mm256_mul_ps(fvec, b_scale_vec);  // 使用b_scale_vec
                fvec = _mm256_min_ps(_mm256_max_ps(fvec, min_val), max_val);
                __m256i res32 = _mm256_cvtps_epi32(fvec);
                __m128i res16 = _mm_packus_epi32(_mm256_extracti128_si256(res32, 0),
                                               _mm256_extracti128_si256(res32, 1));
                __m128i res8 = _mm_packus_epi16(res16, res16);
                _mm_storel_epi64((__m128i*)(b + i), res8);
            }
        }
    }

    // 剩余像素处理（三个通道独立处理）
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

    // 创建临时拷贝（使用对齐内存）
    unsigned char* temp = aligned_alloc(64, plane_size * 3);
    memcpy(temp, image->data, plane_size * 3);
    unsigned char* temp_r = temp;
    unsigned char* temp_g = temp + plane_size;
    unsigned char* temp_b = temp + plane_size * 2;

    // 获取各平面指针
    unsigned char* r = image->data;
    unsigned char* g = r + plane_size;
    unsigned char* b = g + plane_size;

    // 并行处理三个通道
    #pragma omp parallel sections
    {
        // R通道处理
        #pragma omp section
        {
            #pragma omp parallel for collapse(2) schedule(dynamic)
            for (int j = y; j < y + h; j++) {
                for (int i = x; i < x + w; i++) {
                    int total = 0;
                    int count = 0;
                    
                    // 邻域遍历
                    for (int dy = -radius; dy <= radius; dy++) {
                        const int ny = j + dy;
                        if (ny < 0 || ny >= height) continue;
                        
                        const int start_x = MAX(i - radius, 0);
                        const int end_x = MIN(i + radius, width - 1);
                        const int len = end_x - start_x + 1;
                        
                        // SIMD向量化部分
                        __m256i vec_sum = _mm256_setzero_si256();
                        int k = 0;
                        for (; k <= len - 8; k += 8) {
                            // 正确加载8个像素
                            const __m128i pixels = _mm_loadl_epi64(
                                (__m128i*)(temp_r + ny * width + start_x + k));
                            const __m256i pixels32 = _mm256_cvtepu8_epi32(pixels);
                            vec_sum = _mm256_add_epi32(vec_sum, pixels32);
                        }
                        
                        // 计算向量部分总和
                        __m128i sum128 = _mm_add_epi32(
                            _mm256_extracti128_si256(vec_sum, 0),
                            _mm256_extracti128_si256(vec_sum, 1));
                        sum128 = _mm_add_epi32(sum128, _mm_shuffle_epi32(sum128, _MM_SHUFFLE(1,0,3,2)));
                        sum128 = _mm_add_epi32(sum128, _mm_shuffle_epi32(sum128, _MM_SHUFFLE(2,3,0,1)));
                        total += _mm_extract_epi32(sum128, 0);

                        // 标量处理剩余像素
                        for (; k < len; k++) {
                            total += temp_r[ny * width + start_x + k];
                        }
                        count += len;
                    }
                    
                    // 计算平均值
                    if (count > 0) {
                        r[j * width + i] = (total + count/2) / count; // 四舍五入
                    }
                }
            }
        }

        // G通道处理（结构同R通道）
        #pragma omp section
        {
            #pragma omp parallel for collapse(2) schedule(dynamic)
            for (int j = y; j < y + h; j++) {
                for (int i = x; i < x + w; i++) {
                    int total = 0;
                    int count = 0;
                    
                    for (int dy = -radius; dy <= radius; dy++) {
                        const int ny = j + dy;
                        if (ny < 0 || ny >= height) continue;
                        
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

        // B通道处理（结构同R通道）
        #pragma omp section
        {
            #pragma omp parallel for collapse(2) schedule(dynamic)
            for (int j = y; j < y + h; j++) {
                for (int i = x; i < x + w; i++) {
                    int total = 0;
                    int count = 0;
                    
                    for (int dy = -radius; dy <= radius; dy++) {
                        const int ny = j + dy;
                        if (ny < 0 || ny >= height) continue;
                        
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
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    if (strcmp(operation, "add") == 0) {
        adjust_brightness(img1, add_value);
    } else if (strcmp(operation, "average") == 0) {
        if (!inputs[1]) {
            fprintf(stderr, "混合需要两个输入文件\n");
            free_bmp(img1);
            return 1;
        }
        img2 = read_bmp(inputs[1], &error);
        clock_gettime(CLOCK_MONOTONIC, &start);
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
    // sleep(30);
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_spent = (end.tv_sec - start.tv_sec) * 1000.0 
                  + (end.tv_nsec - start.tv_nsec) / 1000000.0;
    printf("\n====== 运行统计 ======\n");
    printf("总耗时: %.3f 毫秒\n", time_spent);

    if (!write_bmp(output, result, &error)) {
        fprintf(stderr, "错误: %s\n", error);
        free_bmp(result);
        return 1;
    }

    free_bmp(result);
    return 0;
}