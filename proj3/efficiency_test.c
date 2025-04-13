#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <immintrin.h>
#include <omp.h>

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
    uint8_t* r;
    uint8_t* g;
    uint8_t* b;
    int width;
    int height;
} PlanarData;

typedef struct {
    BITMAPFILEHEADER file_header;
    BITMAPINFOHEADER info_header;
    int width;
    int height;
    unsigned char *interleaved;  // 交错存储
    unsigned char *backup;       // 数据备份
    PlanarData planar;           // 平面存储
} BMPImage;

// 内存对齐分配
void* aligned_malloc(size_t size) {
    void* ptr;
    #ifdef _WIN32
    ptr = _aligned_malloc(size, 64);
    #else
    posix_memalign(&ptr, 64, size);
    #endif
    return ptr;
}

void aligned_free(void* ptr) {
    #ifdef _WIN32
    _aligned_free(ptr);
    #else
    free(ptr);
    #endif
}

void free_bmp(BMPImage *image) {
    if (image) {
        aligned_free(image->interleaved);
        aligned_free(image->backup);
        aligned_free(image->planar.r);
        aligned_free(image->planar.g);
        aligned_free(image->planar.b);
        free(image);
    }
}

BMPImage *read_bmp(const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (!file) return NULL;

    BITMAPFILEHEADER file_header;
    BITMAPINFOHEADER info_header;
    
    fread(&file_header, sizeof(file_header), 1, file);
    fread(&info_header, sizeof(info_header), 1, file);

    const int width = info_header.biWidth;
    const int height = abs(info_header.biHeight);
    const size_t data_size = info_header.biSizeImage;

    // 读取原始数据
    unsigned char *pixel_data = (unsigned char*)aligned_malloc(data_size);
    fseek(file, file_header.bfOffBits, SEEK_SET);
    fread(pixel_data, 1, data_size, file);
    fclose(file);

    // 创建BMP对象
    BMPImage *image = (BMPImage*)malloc(sizeof(BMPImage));
    image->file_header = file_header;
    image->info_header = info_header;
    image->width = width;
    image->height = height;

    // 初始化存储
    const int pixel_count = width * height;
    const size_t interleaved_size = pixel_count * 3;
    
    // 交错存储
    image->interleaved = (unsigned char*)aligned_malloc(interleaved_size);
    memcpy(image->interleaved, pixel_data, interleaved_size);
    
    // 数据备份
    image->backup = (unsigned char*)aligned_malloc(interleaved_size);
    memcpy(image->backup, pixel_data, interleaved_size);

    // 平面存储
    image->planar.r = (uint8_t*)aligned_malloc(pixel_count);
    image->planar.g = (uint8_t*)aligned_malloc(pixel_count);
    image->planar.b = (uint8_t*)aligned_malloc(pixel_count);
    image->planar.width = width;
    image->planar.height = height;

    // 转换到平面存储
    for (int i = 0; i < pixel_count; ++i) {
        image->planar.b[i] = pixel_data[i*3];
        image->planar.g[i] = pixel_data[i*3+1];
        image->planar.r[i] = pixel_data[i*3+2];
    }

    aligned_free(pixel_data);
    return image;
}

void reset_data(BMPImage *image) {
    const size_t size = image->width * image->height * 3;
    memcpy(image->interleaved, image->backup, size);
    
    const int pixel_count = image->width * image->height;
    for(int i=0; i<pixel_count; i++){
        image->planar.r[i] = image->backup[i*3+2];
        image->planar.g[i] = image->backup[i*3+1];
        image->planar.b[i] = image->backup[i*3];
    }
}

// 原生交错存储版本
void adjust_brightness_interleaved(BMPImage *image, int value) {
    unsigned char *p = image->interleaved;
    const int total = image->width * image->height * 3;
    for (int i = 0; i < total; ++i) {
        int new_val = p[i] + value;
        p[i] = (uint8_t)(new_val < 0 ? 0 : (new_val > 255 ? 255 : new_val));
    }
}

// 原生平面存储版本
void adjust_brightness_planar(PlanarData *planar, int value) {
    const int total = planar->width * planar->height;
    for (int i = 0; i < total; ++i) {
        planar->r[i] = (planar->r[i] + value) < 0 ? 0 : 
                      ((planar->r[i] + value) > 255 ? 255 : planar->r[i] + value);
        planar->g[i] = (planar->g[i] + value) < 0 ? 0 : 
                      ((planar->g[i] + value) > 255 ? 255 : planar->g[i] + value);
        planar->b[i] = (planar->b[i] + value) < 0 ? 0 : 
                      ((planar->b[i] + value) > 255 ? 255 : planar->b[i] + value);
    }
}

// 修正后的SIMD+OpenMP优化的平面存储版本
void adjust_brightness_planar_optimized(PlanarData *planar, int value) {
    const int total = planar->width * planar->height;
    const int abs_value = abs(value);
    const __m256i val_vec = _mm256_set1_epi8((uint8_t)abs_value);

    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < total; i += 32) {
        // 加载32像素的R通道
        __m256i r = _mm256_load_si256((__m256i*)(planar->r + i));
        // 加载32像素的G通道
        __m256i g = _mm256_load_si256((__m256i*)(planar->g + i));
        // 加载32像素的B通道
        __m256i b = _mm256_load_si256((__m256i*)(planar->b + i));

        if (value >= 0) {
            r = _mm256_adds_epu8(r, val_vec);
            g = _mm256_adds_epu8(g, val_vec);
            b = _mm256_adds_epu8(b, val_vec);
        } else {
            r = _mm256_subs_epu8(r, val_vec);
            g = _mm256_subs_epu8(g, val_vec);
            b = _mm256_subs_epu8(b, val_vec);
        }

        _mm256_store_si256((__m256i*)(planar->r + i), r);
        _mm256_store_si256((__m256i*)(planar->g + i), g);
        _mm256_store_si256((__m256i*)(planar->b + i), b);
    }

    // 处理剩余像素（总像素数不是32的倍数时）
    int processed = (total / 32) * 32;
    for (int i = processed; i < total; ++i) {
        planar->r[i] = (planar->r[i] + value) < 0 ? 0 : 
                      ((planar->r[i] + value) > 255 ? 255 : planar->r[i] + value);
        planar->g[i] = (planar->g[i] + value) < 0 ? 0 : 
                      ((planar->g[i] + value) > 255 ? 255 : planar->g[i] + value);
        planar->b[i] = (planar->b[i] + value) < 0 ? 0 : 
                      ((planar->b[i] + value) > 255 ? 255 : planar->b[i] + value);
    }
}

// 修正后的SIMD+OpenMP优化的交错存储版本
void adjust_brightness_interleaved_optimized(BMPImage *image, int value) {
    unsigned char *p = image->interleaved;
    const int total = image->width * image->height * 3;
    const int abs_value = abs(value);
    const __m256i val_vec = _mm256_set1_epi8((uint8_t)abs_value);

    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < total; i += 96) { // 每次处理32像素（96字节）
        // 加载前32字节（10像素 + 2字节）
        __m256i chunk0 = _mm256_load_si256((__m256i*)(p + i));
        // 加载中间32字节
        __m256i chunk1 = _mm256_load_si256((__m256i*)(p + i + 32));
        // 加载最后32字节
        __m256i chunk2 = _mm256_load_si256((__m256i*)(p + i + 64));

        if (value >= 0) {
            chunk0 = _mm256_adds_epu8(chunk0, val_vec);
            chunk1 = _mm256_adds_epu8(chunk1, val_vec);
            chunk2 = _mm256_adds_epu8(chunk2, val_vec);
        } else {
            chunk0 = _mm256_subs_epu8(chunk0, val_vec);
            chunk1 = _mm256_subs_epu8(chunk1, val_vec);
            chunk2 = _mm256_subs_epu8(chunk2, val_vec);
        }

        _mm256_store_si256((__m256i*)(p + i), chunk0);
        _mm256_store_si256((__m256i*)(p + i + 32), chunk1);
        _mm256_store_si256((__m256i*)(p + i + 64), chunk2);
    }

    // 处理剩余字节（总字节数不是96的倍数时）
    int processed = (total / 96) * 96;
    for (int i = processed; i < total; ++i) {
        p[i] = (p[i] + value) < 0 ? 0 : 
              ((p[i] + value) > 255 ? 255 : p[i] + value);
    }
}

void benchmark(BMPImage *image, int value) {
    const int warmup = 10;
    const int runs = 100;
    clock_t start, end;
    double times[4] = {0};
    const char *names[] = {
        "Interleaved (Naive)",
        "Planar (Naive)",
        "Interleaved (SIMD+OpenMP)",
        "Planar (SIMD+OpenMP)"
    };

    // 预热缓存
    printf("Warming up...\n");
    for (int i = 0; i < warmup; ++i) {
        adjust_brightness_interleaved(image, value);
        adjust_brightness_planar(&image->planar, value);
        adjust_brightness_interleaved_optimized(image, value);
        adjust_brightness_planar_optimized(&image->planar, value);
        reset_data(image);
    }

    // 测试原生交错存储
    printf("Testing Interleaved (Naive)...\n");
    start = clock();
    for (int i = 0; i < runs; ++i) {
        adjust_brightness_interleaved(image, value);
        reset_data(image);
    }
    end = clock();
    times[0] = (double)(end - start) * 1000 / CLOCKS_PER_SEC / runs;

    // 测试原生平面存储
    printf("Testing Planar (Naive)...\n");
    start = clock();
    for (int i = 0; i < runs; ++i) {
        adjust_brightness_planar(&image->planar, value);
        reset_data(image);
    }
    end = clock();
    times[1] = (double)(end - start) * 1000 / CLOCKS_PER_SEC / runs;

    // 测试优化交错存储
    printf("Testing Interleaved (Optimized)...\n");
    start = clock();
    for (int i = 0; i < runs; ++i) {
        adjust_brightness_interleaved_optimized(image, value);
        reset_data(image);
    }
    end = clock();
    times[2] = (double)(end - start) * 1000 / CLOCKS_PER_SEC / runs;

    // 测试优化平面存储
    printf("Testing Planar (Optimized)...\n");
    start = clock();
    for (int i = 0; i < runs; ++i) {
        adjust_brightness_planar_optimized(&image->planar, value);
        reset_data(image);
    }
    end = clock();
    times[3] = (double)(end - start) * 1000 / CLOCKS_PER_SEC / runs;

    // 打印结果
    printf("\nPerformance Benchmark (Image: %dx%d)\n", image->width, image->height);
    printf("==============================================\n");
    for (int i = 0; i < 4; ++i) {
        printf("%-25s: %7.4f ms/op\n", names[i], times[i]);
    }
    printf("==============================================\n");
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <input.bmp> <brightness_value>\n", argv[0]);
        return 1;
    }

    // 设置OpenMP线程数（根据CPU核心数调整）
    omp_set_num_threads(4);

    BMPImage *image = read_bmp(argv[1]);
    if (!image) {
        printf("Error reading BMP file\n");
        return 1;
    }

    const int value = atoi(argv[2]);
    printf("Using brightness adjustment: %d\n", value);
    printf("Image dimensions: %d x %d\n", image->width, image->height);
    
    benchmark(image, value);
    
    free_bmp(image);
    return 0;
}