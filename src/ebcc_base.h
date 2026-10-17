#ifndef EBCC_BASE_H
#define EBCC_BASE_H

#include <stddef.h>
#include <stdint.h>

/* Base layer: the field, scaled from [lo, hi] to 16 bits, as one JPEG2000
 * codestream, allocated by an error bound or a target compression ratio. */

/* Decodes any single-component 16-bit codestream, including those of earlier
 * EBCC versions. Fills *field if it is not NULL, else allocates it. Returns 0 on failure. */
int ebcc_base_decode(const uint8_t *codestream, size_t size, float lo, float hi,
        float **field, size_t *height, size_t *width);

/* bound is relative to hi - lo; use_fixed_quantization enables nearest-even
 * coefficient rounding and retains all coding passes for the quantization ablation.
 * Returns an allocated codestream or zero. */
size_t ebcc_base_encode(const float *data, size_t components, size_t height, size_t width,
        float lo, float hi, double bound, uint8_t **codestream, int use_fixed_quantization);
/* base_cr is the target ratio relative to the original float32 input. */
size_t ebcc_base_encode_fixed_cr(const float *data, size_t components, size_t height, size_t width,
        float lo, float hi, float base_cr, uint8_t **codestream);
/* nearest selects the fixed-quantization reconstruction rule from the header. */
int ebcc_base_decode_components(const uint8_t *codestream, size_t size, size_t components,
        float lo, float hi, float **field, size_t *height, size_t *width, int nearest);
#endif /* EBCC_BASE_H */
