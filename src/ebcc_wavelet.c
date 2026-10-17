#include "ebcc_wavelet.h"
#include <limits.h>
#include <string.h>

typedef struct {
    double lifting[4], scale;
    unsigned steps;
} wavelet_t;

/* Both families use the same normalization: low-pass DC gain 1 and
 * high-pass Nyquist gain 2. These are not unit-energy-normalized filters. */
static const wavelet_t wavelets[] = {
    /* EBCC_CDF53: predict/update lifting already has the required gains,
     * so the additional normalization scale is 1. */
    {{-0.5, 0.25, 0, 0}, 1.0, 2},
    /* EBCC_CDF97: predict/update/predict/update lifting reaches the same gains
     * by dividing low coefficients and multiplying high ones by scale. */
    {{-1.586134342059924, -0.052980118572961,
       0.882911075530934, 0.443506852043971}, 1.230174104914001, 4}
};

/* One lifting kernel for both families, precisions, and transform directions.
 * Interleaved samples make reflection at odd/even endpoints identical. */
static void lift(double *x, size_t n, const wavelet_t *wavelet, int inverse, int float_steps) {
    for (unsigned s = 0; s < wavelet->steps; ++s) {
        unsigned step = inverse ? wavelet->steps - 1 - s : s;
        double coefficient = inverse ? -wavelet->lifting[step] : wavelet->lifting[step];
        for (size_t i = (step & 1) ? 0 : 1; i < n; i += 2) {
            size_t left = i ? i - 1 : 1;
            size_t right = i + 1 < n ? i + 1 : i - 1;
            double value = x[i] + coefficient * (x[left] + x[right]);
            x[i] = float_steps ? (double)(float)value : value;
        }
    }
}

static void transform_line(double *x, double *scratch, size_t n,
        ebcc_wavelet_t family, int inverse, int float_steps) {
    if (n < 2) return;
    const wavelet_t *wavelet = &wavelets[family];
    size_t low = (n + 1) / 2;
    if (inverse) {
        for (size_t i = 0; i < n; ++i)
            scratch[i] = (i & 1) ? x[low + i / 2] / wavelet->scale
                                : x[i / 2] * wavelet->scale;
        memcpy(x, scratch, n * sizeof(*x));
        lift(x, n, wavelet, 1, float_steps);
    } else {
        lift(x, n, wavelet, 0, float_steps);
        for (size_t i = 0; i < n; ++i)
            scratch[(i & 1) ? low + i / 2 : i / 2] =
                (i & 1) ? x[i] * wavelet->scale : x[i] / wavelet->scale;
        memcpy(x, scratch, n * sizeof(*x));
    }
}

void ebcc_dwt_1d(double *values, double *scratch, size_t n, unsigned levels,
        ebcc_wavelet_t wavelet, int inverse) {
    size_t sizes[sizeof(size_t) * CHAR_BIT];
    unsigned count = 0;
    for (size_t m = n; m > 1 && (!levels || count < levels); m = (m + 1) / 2)
        sizes[count++] = m;
    for (unsigned i = 0; i < count; ++i)
        transform_line(values, scratch, sizes[inverse ? count - 1 - i : i], wavelet, inverse, 0);
}

static void transform_strided(float *values, double *scratch, size_t n,
        size_t stride, ebcc_wavelet_t wavelet, int inverse) {
    for (size_t i = 0; i < n; ++i) scratch[i] = values[i * stride];
    transform_line(scratch, scratch + n, n, wavelet, inverse, 1);
    for (size_t i = 0; i < n; ++i) values[i * stride] = (float)scratch[i];
}

void ebcc_dwt_2d(float *values, double *scratch, size_t rows, size_t cols,
        ebcc_wavelet_t wavelet, int inverse) {
    if (inverse)
        for (size_t x = 0; x < cols; ++x)
            transform_strided(values + x, scratch, rows, cols, wavelet, 1);
    for (size_t y = 0; y < rows; ++y)
        transform_strided(values + y * cols, scratch, cols, 1, wavelet, inverse);
    if (!inverse)
        for (size_t x = 0; x < cols; ++x)
            transform_strided(values + x, scratch, rows, cols, wavelet, 0);
}
