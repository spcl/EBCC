#include "ebcc_internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "log.h"
#ifdef ENABLE_PERF
#include <sys/prctl.h>
#endif

#include <assert.h>

#define MIN(x, y) ((x)<(y)?(x):(y))
#define TRUE 1
#define FALSE 0

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
#define EBCC_STATIC_ASSERT(cond, msg) _Static_assert(cond, msg)
#else
#define EBCC_STATIC_ASSERT_CAT_(a, b) a##b
#define EBCC_STATIC_ASSERT_CAT(a, b) EBCC_STATIC_ASSERT_CAT_(a, b)
#define EBCC_STATIC_ASSERT(cond, msg) typedef char EBCC_STATIC_ASSERT_CAT(ebcc_static_assert_, __LINE__)[(cond) ? 1 : -1]
#endif

typedef struct {
    uint8_t *buffer;
    size_t size;        // size of the buffer (maximum allowed storage)
    size_t length;      // number of data bytes in the buffer
    size_t offset;
} codec_data_buffer_t;

static void codec_data_buffer_init(codec_data_buffer_t* data) {
    const size_t initial_buffer_size = 1024;
    data->buffer = (uint8_t *) malloc(initial_buffer_size);
    data->size = initial_buffer_size;
    data->length = 0;
    data->offset = 0;
}

static void codec_data_buffer_destroy(codec_data_buffer_t* data) {
    free(data->buffer);
}

typedef struct {
    uint8_t magic[4];
    uint32_t version;
    uint32_t ndims;
    uint32_t reserved;
    uint64_t dims[NDIMS];
    uint64_t chunk_dims[NDIMS];
    uint64_t num_chunks;
    uint64_t chunk_size;
} ebcc_chunking_header_t;

EBCC_STATIC_ASSERT(sizeof(float) == 4, "EBCC serialization requires 32-bit float");
EBCC_STATIC_ASSERT(sizeof(ebcc_header_t) == 48, "EBCC header must be fixed-size");
EBCC_STATIC_ASSERT(sizeof(ebcc_chunking_header_t) == 80, "EBCC chunking header must be fixed-size");

static int checked_mul_size_t(size_t a, size_t b, size_t *out) {
    if (a != 0 && b > SIZE_MAX / a) {
        return FALSE;
    }
    *out = a * b;
    return TRUE;
}

static int product_size_t(const size_t values[NDIMS], size_t *out) {
    size_t product = 1;
    for (size_t i = 0; i < NDIMS; ++i) {
        if (!checked_mul_size_t(product, values[i], &product)) {
            return FALSE;
        }
    }
    *out = product;
    return TRUE;
}

static int size_t_to_uint64(size_t value, uint64_t *out) {
    uint64_t converted = (uint64_t) value;
    if ((size_t) converted != value) {
        return FALSE;
    }
    *out = converted;
    return TRUE;
}

static int uint64_to_size_t(uint64_t value, size_t *out) {
    size_t converted = (size_t) value;
    if ((uint64_t) converted != value) {
        return FALSE;
    }
    *out = converted;
    return TRUE;
}

static int size_t_array_to_uint64_array(const size_t src[NDIMS], uint64_t dst[NDIMS]) {
    for (size_t i = 0; i < NDIMS; ++i) {
        if (!size_t_to_uint64(src[i], &dst[i])) {
            return FALSE;
        }
    }
    return TRUE;
}

static int uint64_array_to_size_t_array(const uint64_t src[NDIMS], size_t dst[NDIMS]) {
    for (size_t i = 0; i < NDIMS; ++i) {
        if (!uint64_to_size_t(src[i], &dst[i])) {
            return FALSE;
        }
    }
    return TRUE;
}

/* Older chunk streams stacked components into rows; retain their geometry
 * validation for decoding only. New encodes never use that representation. */
static int legacy_dims_are_valid(const size_t dims[NDIMS]) {
    size_t height;
    return dims[0] && dims[1] && checked_mul_size_t(dims[0], dims[1], &height) &&
        height >= EBCC_MIN_INTERNAL_IMAGE_DIM && height <= EBCC_MAX_INTERNAL_IMAGE_DIM &&
        dims[2] >= EBCC_MIN_INTERNAL_IMAGE_DIM && dims[2] <= EBCC_MAX_INTERNAL_IMAGE_DIM;
}

static size_t ceil_div_size_t(size_t numerator, size_t denominator) {
    return numerator / denominator + (numerator % denominator != 0);
}

static size_t flat_offset(const size_t indices[NDIMS], const size_t dims[NDIMS]) {
    size_t offset = indices[0];
    for (size_t i = 1; i < NDIMS; ++i) {
        offset = offset * dims[i] + indices[i];
    }
    return offset;
}

static void chunk_origin_from_linear(size_t chunk_linear, const size_t chunk_counts[NDIMS],
        const size_t chunk_dims[NDIMS], size_t origin[NDIMS]) {
    for (size_t dim = NDIMS; dim-- > 0;) {
        size_t chunk_index = chunk_linear % chunk_counts[dim];
        chunk_linear /= chunk_counts[dim];
        origin[dim] = chunk_index * chunk_dims[dim];
    }
}

static int chunks_are_contiguous_slabs(const size_t dims[NDIMS], const size_t chunk_dims[NDIMS]) {
    for (size_t dim = 1; dim < NDIMS; ++dim) {
        if (chunk_dims[dim] != dims[dim]) {
            return FALSE;
        }
    }
    return TRUE;
}

static int chunk_fully_in_bounds(const size_t origin[NDIMS], const size_t dims[NDIMS],
        const size_t chunk_dims[NDIMS]) {
    for (size_t dim = 0; dim < NDIMS; ++dim) {
        if (origin[dim] > dims[dim] || chunk_dims[dim] > dims[dim] - origin[dim]) {
            return FALSE;
        }
    }
    return TRUE;
}

static void copy_chunk_from_data_padded(const float *data, float *chunk_buffer, const size_t dims[NDIMS],
        const size_t origin[NDIMS], const size_t chunk_dims[NDIMS], size_t chunk_size) {
    size_t indices[NDIMS];
    for (size_t linear = 0; linear < chunk_size; ++linear) {
        size_t rem = linear;
        for (size_t dim = NDIMS; dim-- > 0;) {
            size_t index = origin[dim] + (rem % chunk_dims[dim]);
            indices[dim] = MIN(index, dims[dim] - 1);
            rem /= chunk_dims[dim];
        }
        chunk_buffer[linear] = data[flat_offset(indices, dims)];
    }
}

static void copy_chunk_to_data_unpadded(const float *chunk_buffer, float *data, const size_t dims[NDIMS],
        const size_t origin[NDIMS], const size_t chunk_dims[NDIMS], size_t chunk_size) {
    size_t indices[NDIMS];
    for (size_t linear = 0; linear < chunk_size; ++linear) {
        size_t rem = linear;
        int in_bounds = TRUE;
        for (size_t dim = NDIMS; dim-- > 0;) {
            indices[dim] = origin[dim] + (rem % chunk_dims[dim]);
            if (indices[dim] >= dims[dim]) {
                in_bounds = FALSE;
            }
            rem /= chunk_dims[dim];
        }
        if (in_bounds) {
            data[flat_offset(indices, dims)] = chunk_buffer[linear];
        }
    }
}

static int codec_data_buffer_append(codec_data_buffer_t *buffer, const void *src, size_t nbytes) {
    if (nbytes > SIZE_MAX - buffer->length) {
        return FALSE;
    }
    size_t required_size = buffer->length + nbytes;

    size_t new_size = buffer->size;
    while (required_size > new_size) {
        if (new_size > SIZE_MAX / 2) {
            return FALSE;
        }
        new_size *= 2;
    }

    if (new_size > buffer->size) {
        uint8_t *new_buffer = (uint8_t *) realloc(buffer->buffer, new_size);
        if (!new_buffer) {
            return FALSE;
        }
        buffer->buffer = new_buffer;
        buffer->size = new_size;
    }

    memcpy(buffer->buffer + buffer->length, src, nbytes);
    buffer->length += nbytes;
    buffer->offset = buffer->length;
    return TRUE;
}

const char* residual_t_names[] = {
    "FIXED_CR",
    "MAX_ERROR",
    "RELATIVE_ERROR"
};

void print_config(codec_config_t *config) {
    log_info("dimensions:\t(%lu, %lu, %lu)", config->dims[0], config->dims[1], config->dims[2]);
    log_info("chunk dimensions:\t(%lu, %lu, %lu)", config->chunk_dims[0], config->chunk_dims[1], config->chunk_dims[2]);
    log_info("base_cr:\t%f", config->base_cr);
    log_info("residual type:\t%s", residual_t_names[config->residual_compression_type]);
    switch (config->residual_compression_type) {
        case FIXED_CR:
            break;
        case MAX_ERROR:
            log_info("max error:\t%f", config->error);
            break;
        case RELATIVE_ERROR:
            log_info("relative error:\t%f", config->error);
            break;
    }
}

void log_set_level_from_env() {
    const char *env_log_level = getenv("EBCC_LOG_LEVEL");
#ifdef DEBUG
    log_set_level(0); /* Default to TRACE */
    log_info("DEBUG flag enabled in compilation");
#else
    log_set_level(3); /* Default to WARN */
#endif
    if (env_log_level) {
        char *endptr;
        int log_level = strtol(env_log_level, &endptr, 10);
        if (*endptr != '\0') log_warn("Ignore log level: %s, should be in [0, 5]: 0 - TRACE, 1 - DEBUG, 2 - INFO, 3 - WARN, 4 - ERROR, 5 - FATAL", env_log_level);
        else {
            log_set_level(log_level);
            log_info("Set log level to %s", log_level_string(log_level));
        }
    }
}

/* Range of the data; FALSE if it holds a NaN or Inf. */
static int finite_range(const float *data, size_t count, float *minval, float *maxval) {
    *minval = *maxval = data[0];
    for (size_t i = 0; i < count; ++i) {
        if (!isfinite(data[i])) {
            log_error("NaN or Inf found in data at index %zu", i);
            return FALSE;
        }
        if (data[i] < *minval) {
            *minval = data[i];
        }
        if (data[i] > *maxval) {
            *maxval = data[i];
        }
    }
    return TRUE;
}

size_t ebcc_encode(float *data, codec_config_t *config, uint8_t **out_buffer) {
    if (out_buffer) *out_buffer = NULL;
    if (!data || !config || !out_buffer) return 0;
    assert_endianness();
    log_set_level_from_env();

    size_t tot_size = 0;
    if (!ebcc_dims_valid(config->dims) || !product_size_t(config->dims, &tot_size)) {
        log_error("Invalid EBCC dimensions: 1..%d components and spatial sides %d..%d required",
                EBCC_MAX_MCT_COMPONENTS, EBCC_MIN_INTERNAL_IMAGE_DIM, EBCC_MAX_INTERNAL_IMAGE_DIM);
        return 0;
    }
    if (config->residual_compression_type != FIXED_CR && config->residual_compression_type != MAX_ERROR &&
            config->residual_compression_type != RELATIVE_ERROR) {
        log_error("Encoding requires FIXED_CR, MAX_ERROR or RELATIVE_ERROR mode");
        return 0;
    }
    if (config->residual_compression_type == FIXED_CR && (!isfinite(config->base_cr) || config->base_cr <= 0)) {
        log_error("FIXED_CR mode requires a positive finite base_cr");
        return 0;
    }
    print_config(config);

    float minval, maxval;
    if (!finite_range(data, tot_size, &minval, &maxval)) {
        return 0;
    }
    if (minval == maxval) {
        return ebcc_encode_constant(minval, config->dims, out_buffer);
    }

    double range = (double) maxval - minval;
    double error_bound = 0;
    if (config->residual_compression_type != FIXED_CR)
        error_bound = config->residual_compression_type == RELATIVE_ERROR ? config->error * range : config->error;
    if (!isfinite((float) range) || (config->residual_compression_type != FIXED_CR &&
            (!isfinite(error_bound) || error_bound <= 0))) {
        log_error("Encoding requires a positive finite error bound and a finite data range");
        return 0;
    }

#ifdef ENABLE_PERF
    prctl(PR_TASK_PERF_EVENTS_ENABLE);
#endif
    size_t size = config->residual_compression_type == FIXED_CR ?
        ebcc_encode_base_only(data, config->dims, minval, maxval, config->base_cr, out_buffer) :
        ebcc_encode_field(data, config->dims, minval, maxval, error_bound, out_buffer);
#ifdef ENABLE_PERF
    prctl(PR_TASK_PERF_EVENTS_DISABLE);
#endif
    return size;
}

size_t ebcc_encode_chunking(float *data, codec_config_t *config, uint8_t **out_buffer) {
    if (out_buffer) *out_buffer = NULL;
    if (!data || !config || !out_buffer) return 0;
    assert_endianness();
    log_set_level_from_env();

    size_t chunk_dims[NDIMS];
    int chunk_dims_all_zero = TRUE;
    for (size_t i = 0; i < NDIMS; ++i) {
        chunk_dims[i] = config->chunk_dims[i];
        if (chunk_dims[i] != 0) {
            chunk_dims_all_zero = FALSE;
        }
    }
    if (chunk_dims_all_zero) {
        for (size_t i = 0; i < NDIMS; ++i) {
            chunk_dims[i] = config->dims[i];
        }
    }
    if (!ebcc_dims_valid(chunk_dims)) {
        log_fatal("Invalid chunking dimensions: product(chunk_dims[0..%d]) and chunk_dims[%d] must be between %d and %d",
                NDIMS - 2, NDIMS - 1, EBCC_MIN_INTERNAL_IMAGE_DIM, EBCC_MAX_INTERNAL_IMAGE_DIM);
        return 0;
    }

    size_t chunk_counts[NDIMS];
    for (size_t i = 0; i < NDIMS; ++i) {
        if (config->dims[i] == 0 || chunk_dims[i] == 0) {
            log_fatal("Invalid chunking dimensions: dims and chunk_dims must be non-zero");
            return 0;
        }
        chunk_counts[i] = ceil_div_size_t(config->dims[i], chunk_dims[i]);
    }

    size_t chunk_size = 0, num_chunks = 0, total_size = 0, chunk_bytes = 0;
    if (!product_size_t(chunk_dims, &chunk_size) ||
            !product_size_t(chunk_counts, &num_chunks) ||
            !product_size_t(config->dims, &total_size) ||
            !checked_mul_size_t(chunk_size, sizeof(float), &chunk_bytes)) {
        log_fatal("Invalid chunking dimensions: size overflow");
        return 0;
    }
    size_t padded_size = 0;
    if (!checked_mul_size_t(chunk_size, num_chunks, &padded_size)) {
        log_fatal("Invalid chunking dimensions: padded size overflow");
        return 0;
    }
    if (padded_size > total_size && padded_size - total_size > total_size / 10) {
        log_warn("Chunk padding adds %lu values over %lu real values (%.2f%%)",
                padded_size - total_size, total_size,
                ((double) (padded_size - total_size) / (double) total_size) * 100.0);
    }
    int reuse_chunk_data = chunks_are_contiguous_slabs(config->dims, chunk_dims);

    codec_data_buffer_t output;
    codec_data_buffer_init(&output);

    ebcc_chunking_header_t header = {0};
    memcpy(header.magic, EBCC_CHUNKING_HEADER_MAGIC, 4);
    header.version = EBCC_CHUNKING_HEADER_VERSION;
    header.ndims = NDIMS;
    if (!size_t_array_to_uint64_array(config->dims, header.dims) ||
            !size_t_array_to_uint64_array(chunk_dims, header.chunk_dims) ||
            !size_t_to_uint64(num_chunks, &header.num_chunks) ||
            !size_t_to_uint64(chunk_size, &header.chunk_size)) {
        log_fatal("Invalid chunking dimensions: serialized size exceeds uint64_t");
        codec_data_buffer_destroy(&output);
        return 0;
    }

    if (!codec_data_buffer_append(&output, &header, sizeof(header))) {
        log_fatal("Failed to allocate chunked EBCC output header");
        codec_data_buffer_destroy(&output);
        return 0;
    }

    float *chunk_buffer = (float *) malloc(chunk_bytes);
    if (!chunk_buffer) {
        log_fatal("Failed to allocate chunk buffer");
        codec_data_buffer_destroy(&output);
        return 0;
    }

    codec_config_t chunk_config = *config;
    memcpy(chunk_config.dims, chunk_dims, sizeof(chunk_config.dims));
    for (size_t i = 0; i < NDIMS; ++i) {
        chunk_config.chunk_dims[i] = 0;
    }

    for (size_t chunk_linear = 0; chunk_linear < num_chunks; ++chunk_linear) {
        size_t origin[NDIMS];
        chunk_origin_from_linear(chunk_linear, chunk_counts, chunk_dims, origin);

        float *chunk_data = chunk_buffer;
        if (reuse_chunk_data && chunk_fully_in_bounds(origin, config->dims, chunk_dims)) {
            chunk_data = data + flat_offset(origin, config->dims);
        } else {
            copy_chunk_from_data_padded(data, chunk_buffer, config->dims, origin, chunk_dims, chunk_size);
        }

        uint8_t *chunk_stream = NULL;
        size_t chunk_stream_size = ebcc_encode(chunk_data, &chunk_config, &chunk_stream);
        if (chunk_stream_size == 0 || !chunk_stream) {
            log_fatal("Failed to encode chunk %lu", chunk_linear);
            free_buffer(chunk_stream);
            free(chunk_buffer);
            codec_data_buffer_destroy(&output);
            return 0;
        }

        uint64_t chunk_stream_size_u64 = 0;
        if (!size_t_to_uint64(chunk_stream_size, &chunk_stream_size_u64)) {
            log_fatal("Failed to append chunk %lu to chunked EBCC output: chunk size exceeds uint64_t",
                    chunk_linear);
            free_buffer(chunk_stream);
            free(chunk_buffer);
            codec_data_buffer_destroy(&output);
            return 0;
        }
        if (!codec_data_buffer_append(&output, &chunk_stream_size_u64, sizeof(chunk_stream_size_u64)) ||
                !codec_data_buffer_append(&output, chunk_stream, chunk_stream_size)) {
            log_fatal("Failed to append chunk %lu to chunked EBCC output", chunk_linear);
            free_buffer(chunk_stream);
            free(chunk_buffer);
            codec_data_buffer_destroy(&output);
            return 0;
        }
        free_buffer(chunk_stream);
    }

    free(chunk_buffer);
    assert(chunk_size != 0);
    *out_buffer = output.buffer;
    return output.length;
}

size_t ebcc_encode_chunking_compat(float *data, codec_config_t *config, uint8_t **out_buffer) {
    if (out_buffer) *out_buffer = NULL;
    if (!data || !config || !out_buffer) return 0;
    assert_endianness();
    assert(NDIMS == 3);
    log_set_level_from_env();

    codec_config_t compat_config = *config;
    int chunk_dims_all_zero = TRUE;
    for (size_t i = 0; i < NDIMS; ++i) {
        if (compat_config.chunk_dims[i] != 0) {
            chunk_dims_all_zero = FALSE;
            break;
        }
    }

    if (chunk_dims_all_zero) {
        compat_config.chunk_dims[0] = MIN(compat_config.dims[0], EBCC_MAX_MCT_COMPONENTS);
        compat_config.chunk_dims[1] = compat_config.dims[1] > EBCC_MAX_INTERNAL_IMAGE_DIM ?
                1024 : compat_config.dims[1];
        compat_config.chunk_dims[2] = compat_config.dims[2] > EBCC_MAX_INTERNAL_IMAGE_DIM ?
                1024 : compat_config.dims[2];
        log_info("ebcc_encode_chunking_compat chunk dimensions: (%lu, %lu, %lu)",
                compat_config.chunk_dims[0], compat_config.chunk_dims[1], compat_config.chunk_dims[2]);
    }

    if (compat_config.residual_compression_type == RELATIVE_ERROR) {
        /* The relative bound refers to the range of the whole array, not of each chunk. */
        size_t total_size = 0;
        if (!product_size_t(compat_config.dims, &total_size) || total_size == 0) {
            log_fatal("Invalid EBCC dimensions: size overflow or zero-sized data");
            return 0;
        }
        float minval, maxval;
        if (!finite_range(data, total_size, &minval, &maxval)) {
            return 0;
        }
        compat_config.error *= maxval - minval;
        compat_config.residual_compression_type = MAX_ERROR;
    }

    return ebcc_encode_chunking(data, &compat_config, out_buffer);
}

size_t ebcc_decode_chunking(uint8_t *data, size_t data_size, float **out_buffer) {
    if (!data || !out_buffer) return 0;
    assert_endianness();
    log_set_level_from_env();

    if (data_size < sizeof(ebcc_chunking_header_t) ||
            memcmp(data, EBCC_CHUNKING_HEADER_MAGIC, 4) != 0) {
        return ebcc_decode(data, data_size, out_buffer);
    }

    const uint8_t *iter = data;
    const uint8_t *end = data + data_size;
    ebcc_chunking_header_t header;
    memcpy(&header, iter, sizeof(header));
    iter += sizeof(header);

    if (header.version != EBCC_CHUNKING_HEADER_VERSION) {
        log_fatal("Unsupported EBCC chunking header version: %u", header.version);
        return 0;
    }
    if (header.ndims != NDIMS) {
        log_fatal("Unsupported EBCC chunking dimensionality: %u", header.ndims);
        return 0;
    }

    size_t dims[NDIMS];
    size_t chunk_dims[NDIMS];
    size_t header_chunk_size = 0;
    size_t header_num_chunks = 0;
    if (!uint64_array_to_size_t_array(header.dims, dims) ||
            !uint64_array_to_size_t_array(header.chunk_dims, chunk_dims) ||
            !uint64_to_size_t(header.chunk_size, &header_chunk_size) ||
            !uint64_to_size_t(header.num_chunks, &header_num_chunks)) {
        log_fatal("Invalid chunked EBCC data: header field exceeds local SIZE_MAX");
        return 0;
    }
    if (!ebcc_dims_valid(chunk_dims) && !legacy_dims_are_valid(chunk_dims)) {
        log_fatal("Invalid chunked EBCC data: product(chunk_dims[0..%d]) and chunk_dims[%d] must be between %d and %d",
                NDIMS - 2, NDIMS - 1, EBCC_MIN_INTERNAL_IMAGE_DIM, EBCC_MAX_INTERNAL_IMAGE_DIM);
        return 0;
    }

    size_t chunk_counts[NDIMS];
    for (size_t i = 0; i < NDIMS; ++i) {
        if (dims[i] == 0 || chunk_dims[i] == 0) {
            log_fatal("Invalid chunked EBCC data: dims and chunk_dims must be non-zero");
            return 0;
        }
        chunk_counts[i] = ceil_div_size_t(dims[i], chunk_dims[i]);
    }

    size_t expected_chunk_size = 0, expected_num_chunks = 0, total_size = 0, total_bytes = 0;
    if (!product_size_t(chunk_dims, &expected_chunk_size) ||
            !product_size_t(chunk_counts, &expected_num_chunks) ||
            !product_size_t(dims, &total_size) ||
            !checked_mul_size_t(total_size, sizeof(float), &total_bytes)) {
        log_fatal("Invalid chunked EBCC data: size overflow");
        return 0;
    }
    if (header_chunk_size != expected_chunk_size || header_num_chunks != expected_num_chunks) {
        log_fatal("Invalid chunked EBCC data: inconsistent chunk metadata");
        return 0;
    }
    int chunks_are_contiguous = chunks_are_contiguous_slabs(dims, chunk_dims);

    int caller_buffer = *out_buffer != NULL;
    if (!caller_buffer) *out_buffer = (float *) malloc(total_bytes);
    if (!*out_buffer) {
        log_fatal("Failed to allocate chunked EBCC decode output");
        return 0;
    }

    for (size_t chunk_linear = 0; chunk_linear < header_num_chunks; ++chunk_linear) {
        if ((size_t) (end - iter) < sizeof(uint64_t)) {
            log_fatal("Invalid chunked EBCC data: missing chunk size");
            goto fail;
        }

        uint64_t chunk_stream_size_u64 = 0;
        memcpy(&chunk_stream_size_u64, iter, sizeof(chunk_stream_size_u64));
        iter += sizeof(chunk_stream_size_u64);
        size_t chunk_stream_size = 0;
        if (!uint64_to_size_t(chunk_stream_size_u64, &chunk_stream_size)) {
            log_fatal("Invalid chunked EBCC data: chunk payload size exceeds local SIZE_MAX");
            goto fail;
        }
        if (chunk_stream_size > (size_t) (end - iter)) {
            log_fatal("Invalid chunked EBCC data: truncated chunk payload");
            goto fail;
        }

        float *chunk_buffer = NULL;
        size_t decoded_chunk_size = ebcc_decode((uint8_t *) iter, chunk_stream_size, &chunk_buffer);
        if (decoded_chunk_size != header_chunk_size || !chunk_buffer) {
            log_fatal("Invalid chunked EBCC data: decoded chunk %lu has %lu values, expected %lu",
                    chunk_linear, decoded_chunk_size, header_chunk_size);
            free_buffer(chunk_buffer);
            goto fail;
        }

        size_t origin[NDIMS];
        chunk_origin_from_linear(chunk_linear, chunk_counts, chunk_dims, origin);
        if (chunks_are_contiguous && chunk_fully_in_bounds(origin, dims, chunk_dims)) {
            memcpy(*out_buffer + flat_offset(origin, dims), chunk_buffer,
                    header_chunk_size * sizeof(float));
        } else {
            copy_chunk_to_data_unpadded(chunk_buffer, *out_buffer, dims, origin,
                    chunk_dims, header_chunk_size);
        }
        free_buffer(chunk_buffer);
        iter += chunk_stream_size;
    }

    if (iter != end) {
        log_fatal("Invalid chunked EBCC data: trailing payload bytes");
        goto fail;
    }

    return total_size;
fail:
    if (!caller_buffer) {
        free(*out_buffer);
        *out_buffer = NULL;
    }
    return 0;
}

void free_buffer(void *buffer) {
    if (buffer) {
        free(buffer);
    }
}
