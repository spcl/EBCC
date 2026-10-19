/* Version-2 decoding and public stream-version dispatch. */
#include "ebcc_internal.h"
#include "ebcc_base.h"
#include "ebcc_residual.h"
#include "log.h"
#include <math.h>
#include <stdlib.h>

static int residual_valid(const residual_header_t *header, size_t available) {
    return header->size <= available && header->raw_size <= SIZE_MAX &&
        ((!header->size && !header->raw_size && header->step == 0) ||
         (header->size && header->raw_size && isfinite(header->step) && header->step > 0));
}

static size_t decode(const uint8_t *data, size_t size, float **output) {
    if (size < sizeof(ebcc_header_t)) return 0;
    ebcc_header_t header;
    memcpy(&header, data, sizeof(header));
    size_t n = header.reserved;
    float lo, hi;
    memcpy(&lo, &header.minval_bits, sizeof(lo));
    memcpy(&hi, &header.maxval_bits, sizeof(hi));
    if (header.version != EBCC_HEAD_VERSION || (header.flags & ~EBCC_HEADER_FLAGS_MASK) ||
            !n || n > EBCC_MAX_MCT_COMPONENTS || !isfinite(lo) || !isfinite(hi) ||
            header.compressed_size > size - sizeof(header) ||
            header.tail_size != size - sizeof(header) - header.compressed_size) return 0;
    if (header.flags & EBCC_HEADER_FLAG_CONST_FIELD) {
        if (header.flags != EBCC_HEADER_FLAG_CONST_FIELD || lo != hi ||
                header.coeffs_size || header.residual_minval_bits || header.residual_maxval_bits ||
                header.compressed_size || header.tail_size != sizeof(uint64_t)) return 0;
        uint64_t count;
        memcpy(&count, data + sizeof(header), sizeof(count));
        if (count < n * EBCC_MIN_INTERNAL_IMAGE_DIM * EBCC_MIN_INTERNAL_IMAGE_DIM ||
                count > n * EBCC_MAX_INTERNAL_IMAGE_DIM * EBCC_MAX_INTERNAL_IMAGE_DIM ||
                count % n || count > SIZE_MAX / sizeof(float)) return 0;
        if (!*output) *output = malloc((size_t)count * sizeof(float));
        if (!*output) return 0;
        for (size_t i = 0; i < (size_t)count; ++i) (*output)[i] = lo;
        return (size_t)count;
    }
    if (!(hi > lo) || !isfinite(hi - lo) || !header.tail_size) return 0;
    size_t residual_end = sizeof(header) + (size_t)header.compressed_size;
    size_t directory_size = n > 1 ? n * sizeof(residual_header_t) : 0;
    if (directory_size > header.compressed_size ||
            (n > 1 && (header.residual_minval_bits || header.residual_maxval_bits))) return 0;

    size_t offset = sizeof(header) + directory_size;
    uint64_t raw_size = 0;
    ebcc_residual_t residuals[EBCC_MAX_MCT_COMPONENTS] = {{0}};
    for (size_t c = 0; c < n; ++c) {
        residual_header_t rh;
        if (n == 1) {
            rh.raw_size = header.coeffs_size;
            rh.size = header.compressed_size;
            memcpy(&rh.step, &header.residual_minval_bits, sizeof(rh.step));
        } else {
            memcpy(&rh, data + sizeof(header) + c * sizeof(rh), sizeof(rh));
        }
        if (!residual_valid(&rh, residual_end - offset) || rh.raw_size > UINT64_MAX - raw_size) return 0;
        residuals[c] = (ebcc_residual_t){(uint8_t *)data + offset,
                (size_t)rh.size, (size_t)rh.raw_size, rh.step};
        offset += (size_t)rh.size;
        raw_size += rh.raw_size;
    }
    if (offset != residual_end || raw_size != header.coeffs_size) return 0;

    float *field = NULL;
    size_t h = 0, w = 0;
    int ok = ebcc_base_decode_components(data + residual_end, (size_t)header.tail_size,
            n, lo, hi, &field, &h, &w, (header.flags & EBCC_HEADER_FLAG_BASE_FIXED_QUANTIZATION) != 0);
    size_t dims[3] = {n, h, w}, plane = h * w;
    ok = ok && ebcc_dims_valid(dims);
    for (size_t c = 0; ok && c < n; ++c)
        ok = ebcc_residual_decode_mode(&residuals[c], h, w, field + c * plane,
                (header.flags & EBCC_HEADER_FLAG_RESIDUAL_SAMPLE) != 0);
    if (!ok) {
        free(field);
        return 0;
    }
    if (*output) {
        memcpy(*output, field, n * plane * sizeof(float));
        free(field);
    } else {
        *output = field;
    }
    return n * plane;
}

size_t ebcc_decode(uint8_t *data, size_t data_size, float **out_buffer) {
    assert_endianness();
    if (!data || !out_buffer) {
        return 0;
    }

    if (data_size >= 4 && !memcmp(data, EBCC_HEADER_MAGIC, 4)) {
        if (data_size < 5) return 0;
        if (data[4] == EBCC_HEAD_VERSION) return decode(data, data_size, out_buffer);
        if (data[4] != 1) {
            log_error("Unsupported EBCC header version: %u", data[4]);
            return 0;
        }
    }
    return ebcc_decode_ver_0_1(data, data_size, out_buffer);
}
