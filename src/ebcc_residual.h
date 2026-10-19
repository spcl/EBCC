#ifndef EBCC_RESIDUAL_H
#define EBCC_RESIDUAL_H

#include <stddef.h>
#include <stdint.h>

/* Residual layer: a sparse wavelet correction that brings the base
 * reconstruction within a maximum pointwise error of the data.
 *
 * The part of the base error that exceeds a deadzone of 0.8 bound is
 * transformed with one level of floating-point CDF 5/3 lifting, quantized
 * uniformly and stored as (gap, zigzag value) LEB128 records under zstd.
 * Start at s = 8 * (1 - 0.8) * bound. After a failed reconstruction with
 * maximum error E, multiply s by clip((1 - 0.8) * bound / (E - 0.8 * bound),
 * 0.1, 0.8). Accept only a verified float32 reconstruction within the bound,
 * or report failure after 16 trials. */

typedef struct {
    uint8_t *bytes;   /* zstd frame, owned by the caller; NULL when no correction is needed */
    size_t size;      /* bytes in the zstd frame */
    size_t raw_size;  /* bytes of the records inside the frame */
    double step;      /* quantization step in data units; 0 when no correction is needed */
} ebcc_residual_t;

/* Returns 0 if no correction within `bound` could be built. */
int ebcc_residual_encode(const float *data, const float *base, size_t height, size_t width,
        double bound, ebcc_residual_t *residual);

/* Adds the correction to `field`, which holds the base reconstruction. Returns 0 on a corrupt stream. */
int ebcc_residual_decode(const ebcc_residual_t *residual, size_t height, size_t width, float *field);

/* Ablation: original sample residual, including nonzero deadzone restoration.
 * Its payload prefixes the compressed records with an explicit double offset. */
int ebcc_residual_encode_mode(const float *data, const float *base, size_t height, size_t width,
        double bound, ebcc_residual_t *residual, int sample);
int ebcc_residual_decode_mode(const ebcc_residual_t *residual, size_t height, size_t width,
        float *field, int sample);

#endif /* EBCC_RESIDUAL_H */
