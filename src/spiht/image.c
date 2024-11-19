#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <image.h>
#include <assert.h>
//#include "wavelet.h"
#include "spiht.h"
#include <math.h>

void loadPGM(const char *filename, size_t *width, size_t *height, unsigned int *maxVal, void **data) {
    FILE *file = fopen(filename, "rb");
    int elem_size;
    if (!file) {
        perror("Unable to open file");
        exit(EXIT_FAILURE);
    }

    char format[3];
    fscanf(file, "%2s", format);
    if (format[0] != 'P' || format[1] != '5') {
        fprintf(stderr, "Invalid PGM file format\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    fscanf(file, "%lu %lu", width, height);
    fscanf(file, "%u", maxVal);
    fgetc(file); // consume the newline character after maxVal

    assert(*maxVal < 65536u);
    elem_size = (*maxVal < 256u) ? 1 : 2;

    *data = malloc((*width) * (*height) * elem_size);
    if (!*data) {
        perror("Unable to allocate memory");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < (*width) * (*height); ++i) {
        if (elem_size == 1) {
            ((uint8_t *)*data)[i] = fgetc(file);
        } else {
            ((uint16_t *)*data)[i] = fgetc(file) << 8 | fgetc(file);
        }
    }
    fclose(file);
}

void writePGM(const char *filename, size_t width, size_t height, unsigned int maxVal, float *data, float dataMin, float dataMax) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        perror("Unable to open file");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "P5\n%lu %lu\n%u\n", width, height, maxVal);
    int elem_size = (maxVal < 256u) ? 1 : 2;
    for (size_t i = 0; i < width * height; ++i) {
        uint16_t val = (uint16_t) roundf((data[i] - dataMin) / (dataMax - dataMin) * maxVal);
        if (elem_size == 1) {
            fputc(val, file);
        } else {
            fputc(val >> 8, file);
            fputc(val & 0xFF, file);
        }
    }
    fclose(file);
}

void transformToFloatArray(const uint16_t *input, double **output, const double scale, size_t size) {
    *output = (double *)malloc(size * sizeof(double));
    if (!*output) {
        perror("Unable to allocate memory for float array");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < size; ++i) {
        (*output)[i] = scale * (double)input[i];
    }
}

void findMinMax(const uint16_t *array, size_t size, uint16_t *min, uint16_t *max) {
    uint16_t min_val = UINT16_MAX, max_val = 0;
    if (size == 0) {
        return;
    }

    for (size_t i = 0; i < size; ++i) {
        if (array[i] < min_val) {
            min_val = array[i];
        }
        if (array[i] > max_val) {
            max_val = array[i];
        }
    }
    *min = min_val;
    *max = max_val;
}

void findMinMaxf(const float *array, size_t size, float *min, float *max) {
    float min_val = INFINITY, max_val = -INFINITY;
    if (size == 0) {
        return;
    }

    for (size_t i = 0; i < size; ++i) {
        if (array[i] < min_val) {
            min_val = array[i];
        }
        if (array[i] > max_val) {
            max_val = array[i];
        }
    }
    *min = min_val;
    *max = max_val;
}

//TODO: scaling 127 or 255 makes difference
//TODO: test higher num_stages (= 6) -> no difference
//TODO: original encoder + new decoder, original decoder + new encoder

int main() {
    const char *filename = "frame.pgm";
    size_t width, height, coeff_size, wavelet_levels = 6;// 3;
    unsigned int maxVal;
    uint16_t *data, min_c=65534, max_c=65534;
    double *floatData, *reconData, *coeff, scale = 1.0, max_error = 0.0;
    float *floatDataf, *floatDataf_norm, *reconDataf, minValf, maxValf, reconMin, reconMax;
    uint8_t *coeff_buf;

    loadPGM(filename, &width, &height, &maxVal, (void**)&data);
    assert((maxVal >= 256u) && (maxVal < 65536u));
    findMinMax(data, width*height, &min_c, &max_c);
    transformToFloatArray(data, &floatData, scale, width*height);

    floatDataf = (float *)malloc(width * height * sizeof(float));
    for (size_t i = 0; i < width * height; ++i) {
        floatDataf[i] = (float)floatData[i];
    }
    findMinMaxf(floatDataf, width*height, &minValf, &maxValf);
    floatDataf_norm = (float *)malloc(width * height * sizeof(float));
    for (size_t i = 0; i < width * height; ++i) {
        floatDataf_norm[i] = (floatDataf[i] - minValf) / (maxValf - minValf);
    }


    // Use the data array as needed
    // For example, print the first pixel value
    printf("First pixel value: %u\n", data[0]);
    printf("Width: %lu, Height: %lu, MaxVal: %u\n", width, height, maxVal);
    printf("Actual Min: %u, Max: %u\n", min_c, max_c);

    // wavelib_forward_double(floatData, height, width, wavelet_levels, &coeff, &coeff_size);
    spiht_encode(floatDataf_norm, height, width, &coeff_buf, &coeff_size, wavelet_levels);
    //spiht_encode_file("frame.pgm", &coeff_buf, &coeff_size, 8.0);
    reconDataf = (float *)calloc(width * height, sizeof(float));
    spiht_decode(coeff_buf, coeff_size, reconDataf, height, width, coeff_size*8);
    //spiht_decode_file(coeff_buf, coeff_size, "recon_imshrinker.pgm");
    // wavelib_backward_double(reconData, height, width, wavelet_levels, coeff);

    for (size_t i = 0; i < width * height; ++i) {
        if (reconDataf[i] < 0.0) {
            reconDataf[i] = 0.0;
        } else if (reconDataf[i] > 1.0) {
            reconDataf[i] = 1.0;
        }
        reconDataf[i] = reconDataf[i] * (maxValf - minValf) + minValf;
    }

    findMinMaxf(reconDataf, width*height, &reconMin, &reconMax);

    for (size_t i = 0; i < width * height; ++i) {
        float error = fabsf(reconDataf[i] - floatDataf[i]);
        max_error = (max_error > error) ? max_error : error;
    }

    writePGM("recon.pgm", width, height, maxVal, reconDataf, reconMin, reconMax);
    writePGM("frame_f.pgm", width, height, maxVal, floatDataf, 0, maxValf);
    printf("Length of wavelet coefficients: %lu\nLength of image data: %lu\n", coeff_size, height * width);
    printf("Max wavelet transform error: %f, maxValf: %f, minValf: %f\n", max_error, maxValf, minValf);
    printf("MaxVal: %f, Min: %f\n", reconMax, reconMin);

    free(data);
    free(floatDataf);
    free(floatDataf_norm);
    free(reconDataf);
    free(coeff_buf);
    return 0;
}