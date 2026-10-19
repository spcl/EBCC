#ifndef EBCC_WAVELET_H
#define EBCC_WAVELET_H
#include <stddef.h>

typedef enum { EBCC_CDF53, EBCC_CDF97 } ebcc_wavelet_t;

/* Whole-sample symmetric lifting, packed low coefficients before high ones.
 * Both families use low-pass DC gain 1 and high-pass Nyquist gain 2.
 * CDF 5/3 needs scale 1; CDF 9/7 needs an additional lifting normalization.
 * Scratch buffers must not overlap values. Length-one lines are unchanged. */

/* Double-precision 1D transform, recursively applied to the low prefix.
 * levels=0 means full decomposition; scratch holds at least n doubles. */
void ebcc_dwt_1d(double *values, double *scratch, size_t n, unsigned levels,
        ebcc_wavelet_t wavelet, int inverse);

/* One separable spatial level, rows then columns (reversed for inverse).
 * Round each lifting step to float32, preserving residual stream arithmetic.
 * scratch holds at least 2*max(rows, cols) doubles. */
void ebcc_dwt_2d(float *values, double *scratch, size_t rows, size_t cols,
        ebcc_wavelet_t wavelet, int inverse);
#endif
