#include "ebcc_mct.h"
#include "ebcc_wavelet.h"

/* Multilevel CDF 9/7 analysis matrix. Whole-sample symmetric lifting,
 * JPEG2000 unit-DC normalization, and low-prefix recursion handle odd sizes
 * without padding. Unit vectors supply the row-major matrix to opj_set_MCT. */
int ebcc_mct_matrix(size_t n, float *matrix) {
    if (!matrix || n < 2 || n > EBCC_MAX_MCT_COMPONENTS) return 0;
    for (size_t col = 0; col < n; ++col) {
        double x[EBCC_MAX_MCT_COMPONENTS] = {0}, tmp[EBCC_MAX_MCT_COMPONENTS];
        x[col] = 1;
        ebcc_dwt_1d(x, tmp, n, 0, EBCC_CDF97, 0);
        for (size_t row = 0; row < n; ++row) matrix[row*n+col] = (float)x[row];
    }
    return 1;
}
