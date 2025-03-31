#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <immintrin.h>
#include <math.h>

typedef struct __attribute__((aligned(64))) {
    __m256d a;
    __m256d b;
} Block4;

typedef struct {
    int total_pairs;
    int blocks_count;
    Block4* blocks;
    double* remain_a;
    double* remain_b;
    int remain_count;
} OptimizedVector;

double get_timestamp() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

void free_optimized_vector(OptimizedVector* ov) {
    if (ov) {
        free(ov->blocks);
        free(ov->remain_a);
        free(ov->remain_b);
        free(ov);
    }
}

OptimizedVector* preprocess_vectors(const double* a, const double* b, int n) {
    OptimizedVector* ov = malloc(sizeof(OptimizedVector));
    ov->blocks_count = n / 4;
    ov->remain_count = n % 4;

    ov->blocks = aligned_alloc(64, ov->blocks_count * sizeof(Block4));
    ov->remain_a = aligned_alloc(64, 4 * sizeof(double));
    ov->remain_b = aligned_alloc(64, 4 * sizeof(double));

    for (int i = 0; i < ov->blocks_count; ++i) {
        const int base = i * 4;
        ov->blocks[i].a = _mm256_loadu_pd(&a[base]);
        ov->blocks[i].b = _mm256_loadu_pd(&b[base]);
    }

    if (ov->remain_count > 0) {
        memcpy(ov->remain_a, &a[ov->blocks_count*4], ov->remain_count*sizeof(double));
        memcpy(ov->remain_b, &b[ov->blocks_count*4], ov->remain_count*sizeof(double));
    }

    ov->total_pairs = n;
    return ov;
}

double optimized_dot_product(const OptimizedVector* ov) {
    __m256d sum = _mm256_setzero_pd();
    
    for (int i = 0; i < ov->blocks_count; ++i) {
        const __m256d a = ov->blocks[i].a;
        const __m256d b = ov->blocks[i].b;
        sum = _mm256_fmadd_pd(a, b, sum);
    }

    __m128d low = _mm256_extractf128_pd(sum, 0);
    __m128d high = _mm256_extractf128_pd(sum, 1);
    low = _mm_add_pd(low, high);
    double total = _mm_cvtsd_f64(_mm_hadd_pd(low, low));

    for (int i = 0; i < ov->remain_count; ++i) {
        total += ov->remain_a[i] * ov->remain_b[i];
    }
    return total;
}

void read_vectors(double** a, double** b, int* n) {
    char buffer[2][1 << 20];
    if (fgets(buffer[0], sizeof(buffer[0]), stdin) == NULL) {
        fprintf(stderr, "Error: Fail to analyze input\n");
        exit(EXIT_FAILURE);
    }
    buffer[0][strcspn(buffer[0], "\n")] = '\0';
    if (fgets(buffer[1], sizeof(buffer[1]), stdin) == NULL) {
        fprintf(stderr, "Error: Fail to analyze input\n");
        exit(EXIT_FAILURE);
    }
    buffer[1][strcspn(buffer[1], "\n")] = '\0';

    double* vec0 = NULL;
    int capacity = 16, len = 0;
    vec0 = malloc(capacity * sizeof(double));
    char* token = strtok(buffer[0], " ,;");
    while (token) {
        char* end;
        double val = strtod(token, &end);
        if (*end != '\0') {
            fprintf(stderr, "错误：非法字符 '%s'\n", token);
            exit(EXIT_FAILURE);
        }
        
        if (len >= capacity) {
            capacity *= 2;
            vec0 = realloc(vec0, capacity * sizeof(double));
        }
        vec0[len++] = val;
        token = strtok(NULL, " ,;");
    }

    double* vec1 = malloc(len * sizeof(double));
    int idx = 0;
    token = strtok(buffer[1], " ,;");
    while (token) {
        char* end;
        double val = strtod(token, &end);
        if (*end != '\0') {
            fprintf(stderr, "Error: Invalid character '%s'\n", token);
            exit(EXIT_FAILURE);
        }
        
        if (idx >= len) {
            fprintf(stderr, "Error: Vector length not match\n");
            exit(EXIT_FAILURE);
        }
        vec1[idx++] = val;
        token = strtok(NULL, " ,;");
    }

    if (idx != len) {
        fprintf(stderr, "Error: Vector length not match\n");
        exit(EXIT_FAILURE);
    }

    *a = vec0;
    *b = vec1;
    *n = len;
}

int main() {
    double *a, *b;
    int n;
    read_vectors(&a, &b, &n);
    
    OptimizedVector* ov = preprocess_vectors(a, b, n);
    
    const double start = get_timestamp();
    const double result = optimized_dot_product(ov);
    const double elapsed = get_timestamp() - start;
    
    printf("Dot product: %.2f\n", result);
    printf("Execution time: %.9f seconds\n", elapsed);
    
    free(a);
    free(b);
    free_optimized_vector(ov);
    return 0;
}