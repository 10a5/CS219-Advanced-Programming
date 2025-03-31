#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>
#include <immintrin.h>
#include <math.h>

double get_timestamp() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

typedef struct {
    int length;   
    double* data;  
} Vector;

void free_vector(Vector* vec) {
    if (vec) {
        free(vec->data);
        free(vec);
    }
}

Vector* read_vector(const char* prompt) {
    char buffer[5000007]; 
    printf("%s", prompt);
    if (!fgets(buffer, sizeof(buffer), stdin)) return NULL;
    buffer[strcspn(buffer, "\n")] = '\0'; 
    char* token = strtok(buffer, " ,;");
    double* values = NULL;
    int capacity = 16, length = 0;

    values = malloc(capacity * sizeof(double));
    if (!values) return NULL;

    while (token) {
        char* endptr;
        double val = strtod(token, &endptr);
        
        if (endptr == token || *endptr != '\0') {
            fprintf(stderr, "Error: Invalid element '%s'\n", token);
            free(values);
            return NULL;
        }

        if (length >= capacity) {
            capacity *= 2;
            double* new_ptr = realloc(values, capacity * sizeof(double));
            if (!new_ptr) {
                free(values);
                return NULL;
            }
            values = new_ptr;
        }

        values[length++] = val;
        token = strtok(NULL, " ,;"); 
    }

    Vector* vec = malloc(sizeof(Vector));
    if (!vec) {
        free(values);
        return NULL;
    }
    vec->length = length;
    vec->data = values;

    return vec;
}

double simd_dot_product(const double* a, const double* b, int n) {
    __m256d sum_vec = _mm256_setzero_pd();
    int i;
    
    for (i = 0; i <= n - 4; i += 4) {
        __m256d a_vec = _mm256_loadu_pd(a + i);
        __m256d b_vec = _mm256_loadu_pd(b + i);
        sum_vec = _mm256_add_pd(sum_vec, _mm256_mul_pd(a_vec, b_vec));
    }
    
    __m128d low  = _mm256_extractf128_pd(sum_vec, 0);
    __m128d high = _mm256_extractf128_pd(sum_vec, 1);
    low = _mm_add_pd(low, high);
    __m128d result = _mm_hadd_pd(low, low);
    double sum = _mm_cvtsd_f64(result);

    for (; i < n; ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

int main() {
    Vector* vec1 = read_vector("");
    Vector* vec2 = read_vector("");

    if (!vec1 || !vec2) {
        fprintf(stderr, "Error: Fail to analyze the input\n");
        free_vector(vec1);
        free_vector(vec2);
        return 1;
    }

    if (vec1->length != vec2->length) {
        fprintf(stderr, "Error: Vector lenth not match! (%d vs %d)\n", vec1->length, vec2->length);
        free_vector(vec1);
        free_vector(vec2);
        return 2;
    }
    const double start_time = get_timestamp();
    const double result = simd_dot_product(vec1->data, vec2->data, vec1->length);
    const double elapsed = get_timestamp() - start_time;

    printf("Dot product: %.2f\n", result);
    printf("Execution time: %.9f seconds\n", elapsed);

    free_vector(vec1);
    free_vector(vec2);
    return 0;
}   