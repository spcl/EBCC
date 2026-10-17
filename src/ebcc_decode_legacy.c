/* Read-only support for master version-1 and unversioned EBCC streams. */
#include "ebcc_internal.h"
#include "ebcc_base.h"
#include "spiht_re.h"
#include "zstd.h"
#include "log.h"
#include <stdlib.h>
#define TRUE 1
#define FALSE 0

static inline float bits_to_float(uint32_t bits) {
    float value = 0.0f;
    memcpy(&value, &bits, sizeof(value));
    return value;
}

/* Read-only master-format stream (unversioned or version 1). */
typedef struct {
    int const_field;
    float minval, maxval;
    float residual_minval, residual_maxval;  /* versions 0 and 1: range of the SPIHT residual */
    struct { const uint8_t *bytes; size_t size, raw_size; } residual;
    const uint8_t *tail;                     /* codestream, or the sample count of a constant field */
    size_t tail_size;
} stream_t;

static int read_bytes(const uint8_t **iter, const uint8_t *end, void *dst, size_t nbytes) {
    if ((size_t) (end - *iter) < nbytes) {
        return FALSE;
    }
    memcpy(dst, *iter, nbytes);
    *iter += nbytes;
    return TRUE;
}

static int parse_stream(const uint8_t *data, size_t data_size, stream_t *stream) {
    const uint8_t *iter = data, *end = data + data_size;
    uint64_t coeffs_size = 0, compressed_size = 0, tail_size = 0;
    memset(stream, 0, sizeof(*stream));

    if (data_size >= 4 && memcmp(data, EBCC_HEADER_MAGIC, 4) == 0) {
        ebcc_header_t header;
        if (!read_bytes(&iter, end, &header, sizeof(header))) {
            return FALSE;
        }
        stream->const_field = (header.flags & EBCC_HEADER_FLAG_CONST_FIELD) != 0;
        stream->minval = bits_to_float(header.minval_bits);
        stream->maxval = bits_to_float(header.maxval_bits);
        if (header.version != 1) return FALSE;
        stream->residual_minval = bits_to_float(header.residual_minval_bits);
        stream->residual_maxval = bits_to_float(header.residual_maxval_bits);
        coeffs_size = header.coeffs_size;
        compressed_size = header.compressed_size;
        tail_size = header.tail_size;
    } else {
        /* Unversioned layout: native floats and size_t, the tail takes the rest. */
        size_t native_coeffs_size = 0, native_compressed_size = 0;
        if (!read_bytes(&iter, end, &stream->minval, sizeof(float)) ||
                !read_bytes(&iter, end, &stream->maxval, sizeof(float)) ||
                !read_bytes(&iter, end, &native_coeffs_size, sizeof(size_t)) ||
                !read_bytes(&iter, end, &stream->residual_minval, sizeof(float)) ||
                !read_bytes(&iter, end, &stream->residual_maxval, sizeof(float)) ||
                !read_bytes(&iter, end, &native_compressed_size, sizeof(size_t))) {
            return FALSE;
        }
        stream->const_field = stream->minval == stream->maxval;
        coeffs_size = native_coeffs_size;
        compressed_size = native_compressed_size;
        tail_size = (uint64_t) (end - iter) < compressed_size ? 0 : (uint64_t) (end - iter) - compressed_size;
    }

    uint64_t remaining = (uint64_t) (end - iter);
    if (coeffs_size > SIZE_MAX || compressed_size > remaining || tail_size > remaining - compressed_size) {
        return FALSE;
    }
    stream->residual.bytes = iter;
    stream->residual.size = (size_t) compressed_size;
    stream->residual.raw_size = (size_t) coeffs_size;
    stream->tail = iter + compressed_size;
    stream->tail_size = (size_t) tail_size;
    return TRUE;
}

/* Versions 0 and 1: SPIHT-coded residual, scaled to [residual_minval, residual_maxval]. */
static int add_spiht_residual(const stream_t *stream, size_t height, size_t width, float *field) {
    size_t coeffs_size = stream->residual.raw_size, tot_size = height * width;
    uint8_t *coeffs = (uint8_t *) calloc(coeffs_size, sizeof(uint8_t));
    float *residual = (float *) calloc(tot_size, sizeof(float));
    int ok = coeffs && residual &&
            ZSTD_decompress(coeffs, coeffs_size, stream->residual.bytes, stream->residual.size) == coeffs_size;
    if (ok) {
        spiht_decode(coeffs, coeffs_size, residual, height, width, coeffs_size * 8);
        for (size_t i = 0; i < tot_size; ++i) {
            field[i] += residual[i] * (stream->residual_maxval - stream->residual_minval) + stream->residual_minval;
        }
    }
    free(coeffs);
    free(residual);
    return ok;
}

size_t ebcc_decode_ver_0_1(const uint8_t *data, size_t data_size, float **out_buffer) {
    stream_t stream;
    if (!parse_stream(data, data_size, &stream)) {
        log_error("Invalid encoded data: truncated stream");
        return 0;
    }
    int caller_buffer = *out_buffer != NULL, ok = FALSE;
    size_t height = 0, width = 0, tot_size = 0;
    if (stream.const_field) {
        uint64_t samples = 0;
        if (stream.tail_size >= sizeof(samples) && !stream.residual.size) {
            memcpy(&samples, stream.tail, sizeof(samples));
        }
        tot_size = samples > SIZE_MAX / sizeof(float) ? 0 : (size_t) samples;
        if (tot_size && !*out_buffer) {
            *out_buffer = (float *) malloc(tot_size * sizeof(float));
        }
        ok = tot_size && *out_buffer;
        for (size_t i = 0; ok && i < tot_size; ++i) {
            (*out_buffer)[i] = stream.minval;
        }
    } else if (ebcc_base_decode(stream.tail, stream.tail_size, stream.minval, stream.maxval,
            out_buffer, &height, &width)) {
        tot_size = height * width;
        ok = !stream.residual.size || !stream.residual.raw_size ||
                add_spiht_residual(&stream, height, width, *out_buffer);
    }

    if (!ok) {
        log_error("Invalid encoded data: corrupt EBCC stream");
        if (!caller_buffer) {
            free(*out_buffer);
            *out_buffer = NULL;
        }
        return 0;
    }
    log_info("Compression Ratio: %f", (double) (tot_size * sizeof(float)) / data_size);
    return tot_size;
}
