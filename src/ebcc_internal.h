/* Internal stream layouts shared by encoding and compatibility decoding. */
#ifndef EBCC_INTERNAL_H
#define EBCC_INTERNAL_H
#include "ebcc_codec.h"
#include <assert.h>
#include <string.h>
/* The master's 48-byte layout is retained exactly. Version 2 uses reserved
 * for the component count; a single-plane quantization step occupies the two
 * residual range fields. Multi-plane residual metadata lives in the payload. */
typedef struct {
    uint8_t magic[4];
    uint8_t version;
    uint8_t flags;
    uint16_t reserved;
    uint32_t minval_bits;
    uint32_t maxval_bits;
    uint64_t coeffs_size;
    uint32_t residual_minval_bits;
    uint32_t residual_maxval_bits;
    uint64_t compressed_size;
    uint64_t tail_size;
} ebcc_header_t;

typedef struct { uint64_t raw_size, size; double step; } residual_header_t;
typedef char header_layout[(sizeof(ebcc_header_t) == 48) ? 1 : -1];
typedef char residual_layout[(sizeof(residual_header_t) == 24) ? 1 : -1];

static inline void assert_endianness(void) {
    uint32_t value = 0x01020304u;
    uint8_t bytes[sizeof(uint32_t)];
    memcpy(bytes, &value, sizeof(bytes));
    assert(bytes[0] == 0x04);
}

static inline int ebcc_dims_valid(const size_t dims[NDIMS]) {
    return dims[0] >= 1 && dims[0] <= EBCC_MAX_MCT_COMPONENTS &&
        dims[1] >= EBCC_MIN_INTERNAL_IMAGE_DIM && dims[1] <= EBCC_MAX_INTERNAL_IMAGE_DIM &&
        dims[2] >= EBCC_MIN_INTERNAL_IMAGE_DIM && dims[2] <= EBCC_MAX_INTERNAL_IMAGE_DIM;
}

size_t ebcc_encode_field(const float *data, const size_t dims[3], float lo, float hi,
        double bound, uint8_t **output);
size_t ebcc_encode_base_only(const float *data, const size_t dims[3], float lo, float hi,
        float base_cr, uint8_t **output);
size_t ebcc_encode_constant(float value, const size_t dims[3], uint8_t **output);
/* Read-only decoding of master version-1 and unversioned streams. */
size_t ebcc_decode_ver_0_1(const uint8_t *data, size_t data_size, float **out_buffer);
#endif
