#include "ebcc_residual.h"
#include "ebcc_wavelet.h"
#include "ebcc_quantize.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "zstd.h"

#define DEADZONE 0.8            /* fraction of the bound left uncorrected */
#define INITIAL_STEP 8.0        /* optimistic CDF 5/3 estimate: sample error is step / 8 */
#define SAMPLE_INITIAL_STEP 2.0 /* sample-space ablation: no inverse wavelet transform */
#define MAX_TRIALS 16           /* step refinements before giving up */
#define SAMPLE_STEP_SAFETY 0.9 /* retain the sample-space ablation search */
#define MIN_SHRINK 0.1
#define MAX_SHRINK 0.8
#define QUANTIZED_LIMIT 0x1p30  /* keeps zigzag values inside 32 bits */
#define ZSTD_LEVEL 22
#define PLANE_ALIGNMENT 16      /* stream format: coefficient planes are padded to multiples of 16 */
#define VARINT_MAX_BYTES 10
#define RECORD_MAX_BYTES (2 * VARINT_MAX_BYTES)

/* Coefficient plane: the field, edge-extended to the padded size. */
typedef struct {
    size_t height, width, rows, cols;
    float *values;
    double *line;
} plane_t;

static size_t padded(size_t n) {
    return (n + PLANE_ALIGNMENT - 1) / PLANE_ALIGNMENT * PLANE_ALIGNMENT;
}

static void plane_free(plane_t *plane) {
    free(plane->values);
    free(plane->line);
}

static int plane_init(plane_t *plane, size_t height, size_t width, int sample) {
    memset(plane, 0, sizeof(*plane));
    if (!height || !width || height > SIZE_MAX - PLANE_ALIGNMENT || width > SIZE_MAX - PLANE_ALIGNMENT) return 0;
    plane->height = height; plane->width = width;
    plane->rows = sample ? height : padded(height);
    plane->cols = sample ? width : padded(width);
    if (plane->rows > SIZE_MAX / RECORD_MAX_BYTES / plane->cols) return 0;
    plane->values = calloc(plane->rows * plane->cols, sizeof(float));
    plane->line = malloc((plane->rows > plane->cols ? plane->rows : plane->cols) * 2 * sizeof(double));
    if (plane->values && plane->line) return 1;
    plane_free(plane);
    return 0;
}

static double max_error(const float *data, const float *base, const plane_t *correction, size_t height, size_t width) {
    double error = 0;
    for (size_t y = 0; y < height; ++y)
        for (size_t x = 0; x < width; ++x) {
            size_t i = y * width + x;
            /* The decoder's float result is what has to meet the bound. */
            float restored = correction ? (float)((double)base[i] + correction->values[y * correction->cols + x]) : base[i];
            if (!isfinite(restored)) return INFINITY;
            error = fmax(error, fabs((double)data[i] - restored));
        }
    return error;
}

static uint8_t *put_varint(uint8_t *p, uint64_t v) {
    for (; v > 0x7f; v >>= 7) *p++ = (uint8_t)(v | 0x80);
    *p++ = (uint8_t)v;
    return p;
}

static int get_varint(const uint8_t **p, const uint8_t *end, uint64_t *v) {
    *v = 0;
    for (unsigned shift = 0; shift < 7 * VARINT_MAX_BYTES; shift += 7) {
        if (*p == end) return 0;
        uint8_t byte = *(*p)++;
        if (shift == 7 * (VARINT_MAX_BYTES - 1) && byte > 1) return 0;
        *v |= (uint64_t)(byte & 0x7f) << shift;
        if (!(byte & 0x80)) return 1;
    }
    return 0;
}

static double quantized_index(float value, double step, int sample) {
    return sample ? ebcc_nearest(value / step) : round(value / step);
}

/* Quantize the coefficients; on success `restored` holds the correction they decode to. */
static int quantize(const plane_t *coefficients, double step, plane_t *restored, int sample, double offset) {
    size_t count = coefficients->rows * coefficients->cols;
    for (size_t i = 0; i < count; ++i) {
        double q = quantized_index(coefficients->values[i], step, sample);
        if (!isfinite(q) || fabs(q) >= QUANTIZED_LIMIT) return 0;
        restored->values[i] = (float)((int32_t)q * step + (sample && q ? copysign(offset, q) : 0));
    }
    if (!sample) ebcc_dwt_2d(restored->values, restored->line, restored->rows, restored->cols, EBCC_CDF53, 1);
    return 1;
}

/* Nonzero quantized coefficients as (zeros skipped, zigzag value) records under zstd. */
static int pack(const plane_t *coefficients, double step, ebcc_residual_t *residual, int sample, double offset) {
    size_t count = coefficients->rows * coefficients->cols, nonzero = 0, next = 0;
    for (size_t i = 0; i < count; ++i) nonzero += quantized_index(coefficients->values[i], step, sample) != 0;
    if (!nonzero) return 1;
    uint8_t *raw = malloc(nonzero * RECORD_MAX_BYTES), *p = raw;
    if (!raw) return 0;
    for (size_t i = 0; i < count; ++i) {
        int64_t q = (int64_t)quantized_index(coefficients->values[i], step, sample);
        if (!q) continue;
        p = put_varint(p, i - next);
        p = put_varint(p, q < 0 ? (uint64_t)(-2 * q - 1) : (uint64_t)(2 * q));
        next = i + 1;
    }
    size_t capacity = ZSTD_compressBound((size_t)(p - raw));
    size_t prefix = sample ? sizeof(offset) : 0;
    residual->bytes = malloc(prefix + capacity);
    residual->raw_size = (size_t)(p - raw);
    residual->size = residual->bytes ? ZSTD_compress(residual->bytes + prefix, capacity, raw, residual->raw_size, ZSTD_LEVEL) : 0;
    residual->step = step;
    free(raw);
    if (!residual->bytes || ZSTD_isError(residual->size)) return 0;
    if (sample) memcpy(residual->bytes, &offset, sizeof(offset));
    residual->size += prefix;
    return 1;
}

int ebcc_residual_encode_mode(const float *data, const float *base, size_t height, size_t width,
        double bound, ebcc_residual_t *residual, int sample) {
    memset(residual, 0, sizeof(*residual));
    if (!isfinite(bound) || bound <= 0) return 0;
    if (max_error(data, base, NULL, height, width) <= bound) return 1;

    plane_t coefficients, restored;
    int ok = 0;
    if (!plane_init(&coefficients, height, width, sample)) return 0;
    if (!plane_init(&restored, height, width, sample)) goto done;
    /* Soft threshold: only the error beyond the deadzone is corrected. */
    double deadzone = bound * DEADZONE;
    for (size_t y = 0; y < coefficients.rows; ++y) {
        size_t sy = y < height ? y : height - 1;
        for (size_t x = 0; x < coefficients.cols; ++x) {
            size_t sx = x < width ? x : width - 1;
            double delta = (double)data[sy * width + sx] - base[sy * width + sx];
            coefficients.values[y * coefficients.cols + x] = (float)copysign(fmax(0.0, fabs(delta) - deadzone), delta);
        }
    }
    if (!sample) ebcc_dwt_2d(coefficients.values, coefficients.line, coefficients.rows, coefficients.cols, EBCC_CDF53, 0);

    /* Quantization noise may use what the deadzone leaves of the bound;
     * shrink the step until the decoded field verifiably meets it. */
    double margin = bound * (1 - DEADZONE);
    double step = margin * (sample ? SAMPLE_INITIAL_STEP : INITIAL_STEP);
    for (unsigned trial = 0; trial < MAX_TRIALS; ++trial) {
        if (!quantize(&coefficients, step, &restored, sample, deadzone)) goto done;
        double error = max_error(data, base, &restored, height, width);
        if (error <= bound) { ok = pack(&coefficients, step, residual, sample, deadzone); break; }
        /* Here error > bound, so error - deadzone exceeds the margin. */
        double shrink = margin / (error - deadzone);
        if (sample) shrink *= SAMPLE_STEP_SAFETY;
        step *= fmax(MIN_SHRINK, fmin(MAX_SHRINK, shrink));
    }
done:
    plane_free(&coefficients);
    plane_free(&restored);
    if (!ok) { free(residual->bytes); memset(residual, 0, sizeof(*residual)); }
    return ok;
}

int ebcc_residual_decode_mode(const ebcc_residual_t *residual, size_t height, size_t width, float *field, int sample) {
    if (!residual->size && !residual->raw_size) return 1;
    size_t prefix = sample ? sizeof(double) : 0;
    double offset = 0;
    if (residual->size <= prefix) return 0;
    if (sample) {
        memcpy(&offset, residual->bytes, sizeof(offset));
        if (!isfinite(offset) || offset < 0) return 0;
    }
    plane_t correction;
    if (!residual->size || !residual->raw_size || !isfinite(residual->step) || residual->step <= 0 ||
        !plane_init(&correction, height, width, sample)) return 0;
    size_t count = correction.rows * correction.cols, next = 0;
    int ok = 0;
    uint8_t *raw = residual->raw_size <= count * RECORD_MAX_BYTES ? malloc(residual->raw_size) : NULL;
    if (!raw || ZSTD_decompress(raw, residual->raw_size, residual->bytes + prefix, residual->size - prefix) != residual->raw_size) goto done;
    for (const uint8_t *p = raw, *end = raw + residual->raw_size; p < end;) {
        uint64_t gap, z;
        if (!get_varint(&p, end, &gap) || !get_varint(&p, end, &z) ||
            next >= count || gap >= count - next || !z || z > INT32_MAX) goto done;
        int32_t q = (z & 1) ? -(int32_t)(z / 2 + 1) : (int32_t)(z / 2);
        correction.values[next + gap] = (float)(q * residual->step + (sample ? copysign(offset, q) : 0));
        next += gap + 1;
    }
    if (!sample) ebcc_dwt_2d(correction.values, correction.line, correction.rows, correction.cols, EBCC_CDF53, 1);
    ok = 1;
    for (size_t y = 0; y < height; ++y)
        for (size_t x = 0; x < width; ++x) {
            float *value = &field[y * width + x];
            *value = (float)((double)*value + correction.values[y * correction.cols + x]);
            ok &= isfinite(*value) != 0;
        }
done:
    free(raw);
    plane_free(&correction);
    return ok;
}

int ebcc_residual_encode(const float *data, const float *base, size_t height, size_t width,
        double bound, ebcc_residual_t *residual) {
    return ebcc_residual_encode_mode(data, base, height, width, bound, residual, 0);
}

int ebcc_residual_decode(const ebcc_residual_t *residual, size_t height, size_t width, float *field) {
    return ebcc_residual_decode_mode(residual, height, width, field, 0);
}
