/* One encoding pipeline for a spatial frame or a multicomponent volume. */
#include "ebcc_internal.h"
#include "ebcc_base.h"
#include "ebcc_residual.h"
#include <stdlib.h>

/* Internal paper ablations: exact 0/1 override the CMake defaults. */
static int ablation_enabled(const char *name, int fallback) {
    const char *value = getenv(name);
    if (value && !strcmp(value, "1")) return 1;
    if (value && !strcmp(value, "0")) return 0;
    return fallback;
}

/* Constants retain master's value-in-header and sample-count-in-tail layout. */
size_t ebcc_encode_constant(float value, const size_t dims[3], uint8_t **out_buffer) {
    ebcc_header_t header = {0};
    memcpy(header.magic, EBCC_HEADER_MAGIC, 4);
    header.version = EBCC_HEAD_VERSION;
    header.flags = EBCC_HEADER_FLAG_CONST_FIELD;
    header.reserved = (uint16_t)dims[0];
    memcpy(&header.minval_bits, &value, sizeof(value));
    header.maxval_bits = header.minval_bits;
    uint64_t samples = (uint64_t)dims[0] * dims[1] * dims[2];
    header.tail_size = sizeof(samples);
    *out_buffer = malloc(sizeof(header) + sizeof(samples));
    if (!*out_buffer) return 0;
    memcpy(*out_buffer, &header, sizeof(header));
    memcpy(*out_buffer + sizeof(header), &samples, sizeof(samples));
    return sizeof(header) + sizeof(samples);
}

static size_t write_field(const size_t dims[3], float lo, float hi,
        const uint8_t *base, size_t base_size, const ebcc_residual_t *residuals,
        uint8_t **output, int use_fixed_quantization, int sample) {
    size_t n = dims[0], directory_size = n > 1 ? n * sizeof(residual_header_t) : 0;
    size_t residual_size = directory_size;
    uint64_t raw_size = 0;
    for (size_t c = 0; c < n; ++c) {
        if (residuals[c].size > SIZE_MAX - residual_size ||
                residuals[c].raw_size > UINT64_MAX - raw_size) return 0;
        residual_size += residuals[c].size;
        raw_size += residuals[c].raw_size;
    }
    if (residual_size > SIZE_MAX - sizeof(ebcc_header_t) ||
            base_size > SIZE_MAX - sizeof(ebcc_header_t) - residual_size) return 0;
    size_t size = sizeof(ebcc_header_t) + residual_size + base_size;
    uint8_t *out = malloc(size);
    if (!out) return 0;
    ebcc_header_t header = {0};
    memcpy(header.magic, EBCC_HEADER_MAGIC, 4);
    header.version = EBCC_HEAD_VERSION;
    header.flags = (use_fixed_quantization ? EBCC_HEADER_FLAG_BASE_FIXED_QUANTIZATION : 0) |
                   (sample ? EBCC_HEADER_FLAG_RESIDUAL_SAMPLE : 0);
    header.reserved = (uint16_t)n;
    memcpy(&header.minval_bits, &lo, sizeof(lo));
    memcpy(&header.maxval_bits, &hi, sizeof(hi));
    header.coeffs_size = raw_size;
    if (n == 1) memcpy(&header.residual_minval_bits, &residuals[0].step, sizeof(double));
    header.compressed_size = residual_size;
    header.tail_size = base_size;
    memcpy(out, &header, sizeof(header));
    size_t offset = sizeof(header) + directory_size;
    for (size_t c = 0; c < n; ++c) {
        if (n > 1) {
            residual_header_t rh = {residuals[c].raw_size, residuals[c].size, residuals[c].step};
            memcpy(out + sizeof(header) + c * sizeof(rh), &rh, sizeof(rh));
        }
        if (residuals[c].size) memcpy(out + offset, residuals[c].bytes, residuals[c].size);
        offset += residuals[c].size;
    }
    memcpy(out + offset, base, base_size);
    *output = out;
    return size;
}

/* Empty residual records keep the existing single-plane and volume formats. */
size_t ebcc_encode_base_only(const float *data, const size_t dims[3], float lo, float hi,
        float base_cr, uint8_t **output) {
    uint8_t *stream = NULL;
    ebcc_residual_t residuals[EBCC_MAX_MCT_COMPONENTS] = {{0}};
    size_t base_size = ebcc_base_encode_fixed_cr(data, dims[0], dims[1], dims[2],
            lo, hi, base_cr, &stream);
    if (!base_size) return 0;
    size_t size = write_field(dims, lo, hi, stream, base_size, residuals, output, 0, 0);
    free(stream);
    return size;
}

/* Geometry, range and absolute bound are validated once by the public API. */
size_t ebcc_encode_field(const float *data, const size_t dims[3], float lo, float hi,
        double bound, uint8_t **output) {
    int use_fixed_quantization = ablation_enabled("EBCC_INTERNAL_FIXED_QUANTIZATION", EBCC_FIXED_QUANTIZATION);
    int sample = ablation_enabled("EBCC_INTERNAL_SAMPLE_RESIDUAL", EBCC_SAMPLE_RESIDUAL);
    size_t n = dims[0], h = dims[1], w = dims[2], plane = h * w;
    size_t bh = 0, bw = 0, size = 0;
    float *base = NULL;
    uint8_t *stream = NULL;
    ebcc_residual_t residuals[EBCC_MAX_MCT_COMPONENTS] = {{0}};
    size_t base_size = ebcc_base_encode(data, n, h, w, lo, hi,
            bound / ((double)hi - lo), &stream, use_fixed_quantization);
    if (!base_size || !ebcc_base_decode_components(stream, base_size, n, lo, hi,
            &base, &bh, &bw, use_fixed_quantization) || bh != h || bw != w) goto done;
    for (size_t c = 0; c < n; ++c) {
        if (!ebcc_residual_encode_mode(data + c * plane, base + c * plane,
                h, w, bound, &residuals[c], sample)) goto done;
    }
    size = write_field(dims, lo, hi, stream, base_size, residuals, output, use_fixed_quantization, sample);
done:
    free(stream);
    free(base);
    for (size_t c = 0; c < n; ++c) free(residuals[c].bytes);
    return size;
}
