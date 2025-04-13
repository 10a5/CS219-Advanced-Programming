#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

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
    size_t plane_size = image->width * image->height;
    unsigned char *r = image->data;
    unsigned char *g = image->data + plane_size;
    unsigned char *b = image->data + plane_size * 2;

    for (size_t i = 0; i < plane_size; i++) {
        r[i] = fmax(fmin(r[i] + value, 255), 0);
        g[i] = fmax(fmin(g[i] + value, 255), 0);
        b[i] = fmax(fmin(b[i] + value, 255), 0);
    }
}

BMPImage *blend_images(BMPImage *img1, BMPImage *img2, char **error) {
    if (img1->width != img2->width || img1->height != img2->height) {
        *error = "图像尺寸不匹配";
        return NULL;
    }

    BMPImage *result = malloc(sizeof(BMPImage));
    *result = *img1;
    size_t plane_size = img1->width * img1->height;
    result->data = malloc(plane_size * 3);
    if (!result->data) {
        *error = "内存分配失败";
        free(result);
        return NULL;
    }

    unsigned char *r1 = img1->data, *r2 = img2->data, *rr = result->data;
    unsigned char *g1 = img1->data + plane_size, *g2 = img2->data + plane_size, *gr = rr + plane_size;
    unsigned char *b1 = img1->data + plane_size*2, *b2 = img2->data + plane_size*2, *br = gr + plane_size;

    for (size_t i = 0; i < plane_size; i++) {
        rr[i] = (r1[i] + r2[i]) / 2;
        gr[i] = (g1[i] + g2[i]) / 2;
        br[i] = (b1[i] + b2[i]) / 2;
    }

    return result;
}

void adjust_contrast(BMPImage *image, float alpha) {
    size_t plane_size = image->width * image->height;
    unsigned char *r = image->data;
    unsigned char *g = image->data + plane_size;
    unsigned char *b = image->data + plane_size * 2;

    for (size_t i = 0; i < plane_size; i++) {
        r[i] = fmax(fmin(alpha * (r[i] - 128) + 128, 255), 0);
        g[i] = fmax(fmin(alpha * (g[i] - 128) + 128, 255), 0);
        b[i] = fmax(fmin(alpha * (b[i] - 128) + 128, 255), 0);
    }
}

void convert_to_grayscale(BMPImage *image) {
    size_t plane_size = image->width * image->height;
    unsigned char *r = image->data;
    unsigned char *g = image->data + plane_size;
    unsigned char *b = image->data + plane_size * 2;

    for (size_t i = 0; i < plane_size; i++) {
        unsigned char gray = 0.299 * r[i] + 0.587 * g[i] + 0.114 * b[i];
        r[i] = g[i] = b[i] = gray;
    }
}

void invert_colors(BMPImage *image) {
    size_t plane_size = image->width * image->height;
    unsigned char *planes[3] = {
        image->data,
        image->data + plane_size,
        image->data + plane_size * 2
    };

    for (int p = 0; p < 3; p++) {
        for (size_t i = 0; i < plane_size; i++) {
            planes[p][i] = 255 - planes[p][i];
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

    memset(result->data, 0, new_plane_size * 3);

    for (int row = 0; row < h; row++) {
        for (int col = 0; col < w; col++) {
            int src_x = x + col;
            int src_y = y + row;
            
            // 源索引检查
            if (src_x >= 0 && src_x < image->width && 
                src_y >= 0 && src_y < image->height) {
                int src_idx = src_y * image->width + src_x;
                int dst_idx = row * w + col;
                
                dst_r[dst_idx] = src_r[src_idx];
                dst_g[dst_idx] = src_g[src_idx];
                dst_b[dst_idx] = src_b[src_idx];
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

    if (degrees == 90 || degrees == -90) {
        result->width = image->height;
        result->height = image->width;
        result->info_header.biWidth = image->height;
        result->info_header.biHeight = image->width;
    }

    size_t new_plane_size = result->width * result->height;
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

    for (int y = 0; y < image->height; y++) {
        for (int x = 0; x < image->width; x++) {
            int src_idx = y * image->width + x;
            int dst_x, dst_y;

            if (degrees == 90) {
                dst_x = image->height - 1 - y;
                dst_y = x;
            } else if (degrees == -90) {
                dst_x = y;
                dst_y = image->width - 1 - x;
            } else { // 180
                dst_x = image->width - 1 - x;
                dst_y = image->height - 1 - y;
            }

            // 添加边界检查
            if (dst_x >= 0 && dst_x < result->width && dst_y >= 0 && dst_y < result->height) {
                int dst_idx = dst_y * result->width + dst_x;
                dst_r[dst_idx] = src_r[src_idx];
                dst_g[dst_idx] = src_g[src_idx];
                dst_b[dst_idx] = src_b[src_idx];
            }
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
    *result = *image;
    result->width = new_width;
    result->height = new_height;
    result->info_header.biWidth = new_width;
    result->info_header.biHeight = new_height;

    size_t new_plane_size = new_width * new_height;
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

    float x_ratio = (float)image->width / new_width;
    float y_ratio = (float)image->height / new_height;

    for (int y = 0; y < new_height; y++) {
        for (int x = 0; x < new_width; x++) {
            int src_x = (int)(x * x_ratio);
            int src_y = (int)(y * y_ratio);
            src_x = src_x < image->width ? src_x : image->width - 1;
            src_y = src_y < image->height ? src_y : image->height - 1;

            int src_idx = src_y * image->width + src_x;
            int dst_idx = y * new_width + x;
            dst_r[dst_idx] = src_r[src_idx];
            dst_g[dst_idx] = src_g[src_idx];
            dst_b[dst_idx] = src_b[src_idx];
        }
    }

    return result;
}

// ============== 新增功能 ==============
void apply_vignette(BMPImage *image, float intensity) {
    size_t plane_size = image->width * image->height;
    unsigned char *r = image->data;
    unsigned char *g = image->data + plane_size;
    unsigned char *b = image->data + plane_size * 2;

    const float center_x = (image->width - 1) / 2.0f;
    const float center_y = (image->height - 1) / 2.0f;
    const float max_radius = sqrtf(center_x*center_x + center_y*center_y);

    for (int y = 0; y < image->height; y++) {
        for (int x = 0; x < image->width; x++) {
            const float dx = x - center_x;
            const float dy = y - center_y;
            const float distance = sqrtf(dx*dx + dy*dy) / max_radius;
            const float factor = 1.0f - intensity * distance * distance;

            const int idx = y * image->width + x;
            r[idx] = (unsigned char)fmaxf(0, fminf(r[idx] * factor, 255));
            g[idx] = (unsigned char)fmaxf(0, fminf(g[idx] * factor, 255));
            b[idx] = (unsigned char)fmaxf(0, fminf(b[idx] * factor, 255));
        }
    }
}

void adjust_rgb_channels(BMPImage *image, float r_scale, float g_scale, float b_scale) {
    size_t plane_size = image->width * image->height;
    unsigned char *r = image->data;
    unsigned char *g = image->data + plane_size;
    unsigned char *b = image->data + plane_size * 2;

    for (size_t i = 0; i < plane_size; i++) {
        r[i] = fmin(r[i] * r_scale, 255);
        g[i] = fmin(g[i] * g_scale, 255);
        b[i] = fmin(b[i] * b_scale, 255);
    }
}

void partial_blur(BMPImage *image, int x, int y, int w, int h, int radius) {
    size_t plane_size = image->width * image->height;
    unsigned char *temp = malloc(plane_size * 3);
    memcpy(temp, image->data, plane_size * 3);
    unsigned char *temp_r = temp;
    unsigned char *temp_g = temp + plane_size;
    unsigned char *temp_b = temp + plane_size * 2;

    unsigned char *r = image->data;
    unsigned char *g = image->data + plane_size;
    unsigned char *b = image->data + plane_size * 2;

    for (int j = y; j < y + h; j++) {
        for (int i = x; i < x + w; i++) {
            int sum_r = 0, sum_g = 0, sum_b = 0, count = 0;
            
            for (int dy = -radius; dy <= radius; dy++) {
                for (int dx = -radius; dx <= radius; dx++) {
                    int nx = i + dx;
                    int ny = j + dy;
                    if (nx >= 0 && nx < image->width && ny >= 0 && ny < image->height) {
                        int idx = ny * image->width + nx;
                        sum_r += temp_r[idx];
                        sum_g += temp_g[idx];
                        sum_b += temp_b[idx];
                        count++;
                    }
                }
            }
            
            if (count > 0) {
                int idx = j * image->width + i;
                r[idx] = sum_r / count;
                g[idx] = sum_g / count;
                b[idx] = sum_b / count;
            }
        }
    }
    free(temp);
}

// 主函数保持不变...
// [注意：主函数部分与原始代码相同，此处省略以节省空间，实际使用时需保留完整主函数]

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
    // clock_t end_time = clock();
    // double time_spent = (double)(end_time - start_time) * 1000.0 / CLOCKS_PER_SEC;
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