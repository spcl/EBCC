#include "wavelib.h"
#include <math.h>

void wavelib_forward(float *image, size_t height, size_t width, size_t levels, double **out, size_t *out_size) {
    size_t size = height * width;
    wave_object obj = wave_init("bior3.3");
    wt2_object wt = wt2_init(obj, "dwt", height, width, levels);
    setDWT2Extension(wt, "per");

    double *wave_buffer = (double *) malloc(size * sizeof(double));
    for (size_t i = 0; i < size; ++i) {
        wave_buffer[i] = image[i];
    }

    *out = dwt2(wt, wave_buffer);
    *out_size = wt->outlength;

    free(wave_buffer);
}

void wavelib_backward(float *image, size_t height, size_t width, size_t levels, double *coeffs) {
    size_t size = height * width;
    wave_object obj = wave_init("bior3.3");
    wt2_object wt = wt2_init(obj, "dwt", height, width, levels);
    setDWT2Extension(wt, "per");

    double *wave_buffer = (double *) malloc(size * sizeof(double));
    free(dwt2(wt, wave_buffer)); // sadly needed
    idwt2(wt, coeffs, wave_buffer);
    for (size_t i = 0; i < size; ++i) {
        image[i] = wave_buffer[i];
    }

    free(wave_buffer);
}

void wavelib_transform(float *data, float *image, size_t height, size_t width, size_t levels, int decode) {
        size_t size = height * width;
        wave_object obj = wave_init("bior3.3");
        wt2_object wt = wt2_init(obj, "dwt", height, width, levels);
        setDWT2Extension(wt, "per");

        double *wave_buffer = (double *) malloc(size * sizeof(double));
        for (size_t i = 0; i < size; ++i) {
            wave_buffer[i] = image[i];
        }

        double *wavecoeffs;

        if (decode) {
            double *fake_wave = (double *) malloc(size * sizeof(double));
            wavecoeffs = dwt2(wt, fake_wave);
            free(fake_wave);
            double *wave_out = (double *) malloc(size * sizeof(double));
            for (size_t i = 0; i < size; ++i) {
                wave_buffer[i] = data[i];
            }
            idwt2(wt, wave_buffer, wave_out);
            for (size_t i = 0; i < size; ++i) {
                image[i] = wave_out[i];
            }
            free(wave_out);
        } else {
            wavecoeffs = dwt2(wt, wave_buffer);
            printf("%d vs %d\n", wt->outlength, (int) size);
            for (size_t i = 0; i < size; ++i) {
                data[i] = wavecoeffs[i];
            }
        }


        free(wavecoeffs);
        free(wave_buffer);
}