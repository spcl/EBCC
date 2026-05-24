#include "ebcc_codec.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "log.h"
#ifdef ENABLE_PERF
#include <sys/prctl.h>
#endif

#include "openjpeg.h"
#include "zstd.h"
#include <assert.h>
#ifdef __linux__
#include <malloc.h>
#else
#include <stdlib.h>
#endif
#include "spiht_re.h"


#define MIN(x, y) ((x)<(y)?(x):(y))
#define CEIL(x,y) (((x) + (y) - 1) / (y))
#define TRUE 1
#define FALSE 0

#define WAVELET_LEVELS 3

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

void j2k_decode_internal(float **data, size_t *height, size_t *width, float minval, float maxval,
        codec_data_buffer_t *codec_data_buffer);

void codec_data_buffer_init(codec_data_buffer_t* data) {
    const size_t initial_buffer_size = 1024;
    data->buffer = (uint8_t *) malloc(initial_buffer_size);
    data->size = initial_buffer_size;
    data->length = 0;
    data->offset = 0;
}

void codec_data_buffer_rewind(codec_data_buffer_t* data) {
    data->offset = 0;
}

void codec_data_buffer_clear(codec_data_buffer_t* data) {
    data->length = 0;
    data->offset = 0;
}

void codec_data_buffer_destroy(codec_data_buffer_t* data) {
    free(data->buffer);
}

void codec_data_buffer_dump(codec_data_buffer_t* data, const char *file_name) {
    FILE *file = fopen(file_name, "wb");
    fwrite(data->buffer, sizeof(char), data->length, file);
    fclose(file);
}

OPJ_SIZE_T write_to_buffer_stream(void *input_buffer, OPJ_SIZE_T len, codec_data_buffer_t *stream_data) {
    size_t new_size = stream_data->size;
    while (stream_data->length + len > new_size) {
        new_size *= 2;
    }

    if (new_size > stream_data->size) {
        stream_data->buffer = (uint8_t *) realloc(stream_data->buffer, new_size);
        stream_data->size = new_size;
    }

    memcpy(stream_data->buffer + stream_data->offset, input_buffer, len);
    stream_data->offset += len;
    stream_data->length += len;

    return len;
}

OPJ_SIZE_T read_from_buffer_stream(void *output_buffer, OPJ_SIZE_T len, codec_data_buffer_t *stream_data) {
    if (stream_data->offset >= stream_data->length) {
        return -1;
    }

    size_t n_bytes_to_read = MIN(len, stream_data->length - stream_data->offset);
    memcpy(output_buffer, stream_data->buffer + stream_data->offset, n_bytes_to_read);
    stream_data->offset += n_bytes_to_read;

    return n_bytes_to_read;
}

void j2k_encode_internal(void *data, size_t *image_dims, size_t *tile_dims, float base_cr,
        codec_data_buffer_t *codec_data_buffer) {
    size_t n_tiles = image_dims[0] / tile_dims[0];
    size_t tile_size = tile_dims[0] * tile_dims[1];

    opj_cparameters_t parameters;
    opj_set_default_encoder_parameters(&parameters);

    // Set image parameters
    parameters.tcp_numlayers = 1;
    parameters.cp_disto_alloc = 1;
    parameters.tcp_rates[0] = base_cr / 2;
    parameters.irreversible = 1;
    parameters.cp_tx0 = 0;
    parameters.cp_ty0 = 0;

    if (n_tiles > 1) {
        parameters.tile_size_on = OPJ_TRUE;
        parameters.cp_tdx = tile_dims[1];
        parameters.cp_tdy = tile_dims[0];
    }

    opj_image_cmptparm_t cmptparm;
    memset(&cmptparm, 0, sizeof(opj_image_cmptparm_t));

    cmptparm.dx = 1;
    cmptparm.dy = 1;
    cmptparm.w = image_dims[1];
    cmptparm.h = image_dims[0];
    cmptparm.prec = 16;
    cmptparm.sgnd = 0;
    cmptparm.x0 = 0;
    cmptparm.y0 = 0;

    opj_image_t *image;
    if (n_tiles == 1) {
        image = opj_image_create(1, &cmptparm, OPJ_CLRSPC_GRAY);
        for (size_t i = 0; i < image_dims[0] * image_dims[1]; ++i) {
            image->comps[0].data[i] = ((uint16_t *) data)[i];
        }
    } else {
        image = opj_image_tile_create(1, &cmptparm, OPJ_CLRSPC_GRAY);
    }
    image->x0 = 0;
    image->y0 = 0;
    image->x1 = image_dims[1];
    image->y1 = image_dims[0];

    opj_codec_t *codec = opj_create_compress(OPJ_CODEC_J2K);

    // Initialize the compressor
    opj_setup_encoder(codec, &parameters, image);

    opj_stream_t *stream = opj_stream_default_create(OPJ_FALSE);

    opj_stream_set_user_data(stream, codec_data_buffer, NULL);
    opj_stream_set_user_data_length(stream, 0);
    opj_stream_set_write_function(stream, (opj_stream_write_fn) write_to_buffer_stream);

    // Compress the image
    opj_start_compress(codec, image, stream);

    if (n_tiles > 1) {
        for (OPJ_UINT32 i = 0; i < n_tiles; ++i) {
            opj_write_tile(codec, i, ((OPJ_BYTE *) data) + i * (tile_size * sizeof(uint16_t)),
                    tile_size * sizeof(uint16_t), stream);
        }
    } else {
        opj_encode(codec, stream);
    }

    opj_end_compress(codec, stream);
    opj_stream_destroy(stream);
    opj_image_destroy(image);
    opj_destroy_codec(codec);
}

typedef uint8_t coo_v_t;

typedef struct {
    uint16_t c;
    coo_v_t v;
} coo_t;


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

static uint32_t float_to_bits(float value) {
    uint32_t bits = 0;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

static float bits_to_float(uint32_t bits) {
    float value = 0.0f;
    memcpy(&value, &bits, sizeof(value));
    return value;
}

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

static int dims_are_valid(const size_t dims[NDIMS]) {
    size_t image_height = 1;
    for (size_t i = 0; i < NDIMS - 1; ++i) {
        if (dims[i] == 0 || !checked_mul_size_t(image_height, dims[i], &image_height)) {
            return FALSE;
        }
    }
    return image_height >= EBCC_MIN_INTERNAL_IMAGE_DIM &&
            image_height <= EBCC_MAX_INTERNAL_IMAGE_DIM &&
            dims[NDIMS - 1] >= EBCC_MIN_INTERNAL_IMAGE_DIM &&
            dims[NDIMS - 1] <= EBCC_MAX_INTERNAL_IMAGE_DIM;
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

static void assert_endianness(void) {
    uint32_t value = 0x01020304u;
    uint8_t bytes[sizeof(uint32_t)];
    memcpy(bytes, &value, sizeof(bytes));
    assert(bytes[0] == 0x04);
}

const char* residual_t_names[] = {
    "NONE",
    "MAX_ERROR",
    "RELATIVE_ERROR"
};

void print_config(codec_config_t *config) {
    log_info("dimensions:\t(%lu, %lu, %lu)", config->dims[0], config->dims[1], config->dims[2]);
    log_info("chunk dimensions:\t(%lu, %lu, %lu)", config->chunk_dims[0], config->chunk_dims[1], config->chunk_dims[2]);
    log_info("base_cr:\t%f", config->base_cr);
    log_info("residual type:\t%s", residual_t_names[config->residual_compression_type]);
    switch (config->residual_compression_type) {
        case NONE:
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

float get_data_range(const float *data, const size_t tot_size) {
    float max_data = data[0], min_data = data[0];
    for (size_t i = 1; i < tot_size; ++i) {
        if (data[i] > max_data) {
            max_data = data[i];
        }
        if (data[i] < min_data) {
            min_data = data[i];
        }
    }
    return max_data - min_data;
}

/*Value-Range Relative error*/
float get_max_relative_error(const float *data, const float *decoded, const float *residual, const size_t tot_size, const float data_range) {
    float cur_max_error = 0;
    assert(tot_size > 0);
    for (size_t i = 0; i < tot_size; ++i) {
        float residual_value = residual ? residual[i] : 0;
        float cur_error = fabsf(data[i] - decoded[i] - residual_value) / data_range;
        if (cur_error > cur_max_error) {
            cur_max_error = cur_error;
        }
    }
    return cur_max_error;
}

float get_max_error(const float *data, const float *decoded, const float *residual, const size_t tot_size) {
    float cur_max_error = 0;
    for (size_t i = 0; i < tot_size; ++i) {
        float residual_value = residual ? residual[i] : 0;
        float cur_error = fabsf(data[i] - (decoded[i] + residual_value));
        /* this is pointwise relative error
        if (error_type == RELATIVE_ERROR) {
            cur_error /= fabsf(data[i]);
        }
        */
        if (cur_error > cur_max_error) {
            cur_max_error = cur_error;
        }
    }
    return cur_max_error;
}

double get_mean_error(const float *data, const float *decoded, const float *residual, const size_t tot_size) {
    double sum_error = 0;
    for (size_t i = 0; i < tot_size; ++i) {
        float residual_value = residual ? residual[i] : 0;
        sum_error += data[i] - (decoded[i] + residual_value);
    }
    return (sum_error / tot_size);
}

double get_error_target_quantile(const float *data, const float *decoded, const float *residual, const size_t tot_size, const float error_target) {
    size_t n = 0;
    for (size_t i = 0; i < tot_size; ++i) {
        float residual_value = residual ? residual[i] : 0;
        float cur_error = fabsf(data[i] - (decoded[i] + residual_value));
        if (cur_error > error_target) {
            n++;
        }
    }
    return 1. -  ((double) n / tot_size);
}

void findMinMaxf(const float *array, size_t size, float *min, float *max) {
    float min_val = INFINITY, max_val = -INFINITY;
    if (size == 0) {
        return;
    }
    min_val = array[0];
    max_val = array[0];

    for (size_t i = 0; i < size; ++i) {
        if (array[i] < min_val) {
            min_val = array[i];
        }
        if (array[i] > max_val) {
            max_val = array[i];
        }
    }
    *min = min_val;
    *max = max_val;
}

double emulate_j2k_compression(uint16_t *scaled_data, size_t *image_dims, size_t *tile_dims, float current_cr, 
                             codec_data_buffer_t *codec_data_buffer, float **decoded, float minval, float maxval, 
                             float *data, size_t tot_size, float error_target) {
    codec_data_buffer_clear(codec_data_buffer);
    j2k_encode_internal(scaled_data, image_dims, tile_dims, current_cr, codec_data_buffer);
    codec_data_buffer_rewind(codec_data_buffer);
    j2k_decode_internal(decoded, NULL, NULL, minval, maxval, codec_data_buffer);
    return get_error_target_quantile(data, *decoded, NULL, tot_size, error_target);
}

float error_bound_j2k_compression(uint16_t *scaled_data, size_t *image_dims, size_t *tile_dims, float current_cr, 
                             codec_data_buffer_t *codec_data_buffer, float **decoded, float minval, float maxval, 
                             float *data, size_t tot_size, float error_target, double base_quantile_target) {
    float cr_lo = current_cr;
    float cr_hi = current_cr;
    double error_target_quantile = get_error_target_quantile(data, *decoded, NULL, tot_size, error_target);
    double error_target_quantile_prev = error_target_quantile;
    double eps = 1e-8;

    log_trace("current_cr: %f, 1-error_target_quantile: %.1e, jp2_length: %lu", current_cr, 1-error_target_quantile, codec_data_buffer->length);

    /* TODO: log down best feasible cr!, best feasible error*/
    /* TODO: take error target quantile from env*/
    /* TODO: log according to env, with log.c */
    while ((error_target_quantile < base_quantile_target) && (cr_lo >= 1./2)) {
        cr_lo /= 2;
        error_target_quantile = emulate_j2k_compression(scaled_data, image_dims, tile_dims, cr_lo, codec_data_buffer, decoded, minval, maxval, data, tot_size, error_target);
        log_trace("cr_lo: %f, 1-error_target_quantile: %.1e, jp2_length: %lu", cr_lo, 1-error_target_quantile, codec_data_buffer->length);
    }
    error_target_quantile = error_target_quantile_prev;
    while ((error_target_quantile >= base_quantile_target) && (cr_hi <= 1000)) {
        cr_hi *= 2;
        error_target_quantile = emulate_j2k_compression(scaled_data, image_dims, tile_dims, cr_hi, codec_data_buffer, decoded, minval, maxval, data, tot_size, error_target);
        log_trace("cr_hi: %f, 1-error_target_quantile: %.1e, jp2_length: %lu", cr_hi, 1-error_target_quantile, codec_data_buffer->length);
    }

    if (error_target_quantile >= base_quantile_target) {
        log_trace("cr_hi: %f exceeds 1000, while base compression still feasible, stay at this level", cr_hi);
        return cr_hi;
    }

    error_target_quantile = error_target_quantile_prev;

    assert(cr_lo <= cr_hi);
    while ((fabs(error_target_quantile - base_quantile_target) > eps || error_target_quantile == 1.0) && cr_hi - cr_lo > 1.) {
        current_cr = (cr_lo + cr_hi) / 2;
        error_target_quantile = emulate_j2k_compression(scaled_data, image_dims, tile_dims, current_cr, codec_data_buffer, decoded, minval, maxval, data, tot_size, error_target);
        log_trace("current_cr: %f, cr_lo: %f, cr_hi: %f, 1-error_target_quantile: %.1e, jp2_length: %lu", current_cr, cr_lo, cr_hi, 1-error_target_quantile, codec_data_buffer->length);
        if (error_target_quantile < base_quantile_target) {
            cr_hi = current_cr;
        } else {
            cr_lo = current_cr;
        }
    }
    /* use cr_lo as the best feasible cr */
    error_target_quantile = emulate_j2k_compression(scaled_data, image_dims, tile_dims, cr_lo, codec_data_buffer, decoded, minval, maxval, data, tot_size, error_target);
    if (error_target_quantile < base_quantile_target) {
        log_warn("Could not reach error target quantile of (1-%.2e) (1-%.2e instead).", 1-base_quantile_target, 1-error_target_quantile);
    }
    return cr_lo;

}

void check_nan_inf(const float *data, const size_t tot_size) {
    for (size_t i = 0; i < tot_size; ++i) {
        if (isnan(data[i]) || isinf(data[i])) {
            log_fatal("NaN or Inf found in data at index %lu", i);
            exit(1);
        }
    }
}

size_t ebcc_encode(float *data, codec_config_t *config, uint8_t **out_buffer) {
#ifdef ENABLE_PERF
    prctl(PR_TASK_PERF_EVENTS_ENABLE);
#endif
    assert_endianness();

    if (!dims_are_valid(config->dims)) {
        log_fatal("Invalid EBCC dimensions: product(dims[0..%d]) and dims[%d] must be between %d and %d",
                NDIMS - 2, NDIMS - 1, EBCC_MIN_INTERNAL_IMAGE_DIM, EBCC_MAX_INTERNAL_IMAGE_DIM);
        return 0;
    }

    int pure_base_codec_required = FALSE, pure_base_codec_done = FALSE, pure_base_codec_disabled = FALSE, pure_base_codec_consistency_disabled = FALSE, mean_error_adjustment_disabled = FALSE;
    size_t compressed_size = 0, base_codec_buffer_length = 0;
    uint8_t *compressed_coefficients = NULL;
    uint8_t *coeffs_buf = NULL;
    uint8_t *base_codec_buffer = NULL; 
    float residual_maxval = 0., residual_minval = 0., error_target = -1, current_cr = -1;
    size_t coeffs_size = 0, coeffs_size_orig = 0, coeffs_trunc_bits = 0; /*coeffs_size: #bytes*/
    double trunc_hi, trunc_lo, best_feasible_trunc;
    double eps = 1e-8, base_error_quantile=1e-6;
    double cur_mean_error = 0.0;

    // Load base_error_quantile from env var EBCC_INIT_BASE_ERROR_QUANTILE, default to 1e-6 if not set

    log_set_level_from_env();
    
    const char *env_base_error_quantile = getenv("EBCC_INIT_BASE_ERROR_QUANTILE");
    const char *env_disable_pure_base_codec_fallback = getenv("EBCC_DISABLE_PURE_BASE_COMPRESSION_FALLBACK");
    const char *env_disable_pure_base_codec_consistency = getenv("EBCC_DISABLE_PURE_BASE_COMPRESSION_FALLBACK_CONSISTENCY");
    const char *env_disable_mean_adjustment = getenv("EBCC_DISABLE_MEAN_ADJUSTMENT");
    if (env_base_error_quantile) {
        base_error_quantile = strtod(env_base_error_quantile, NULL);
    }
    if (env_disable_pure_base_codec_fallback) {
        pure_base_codec_disabled = TRUE;
    }
    if (env_disable_pure_base_codec_consistency) {
        pure_base_codec_consistency_disabled = TRUE;
    }
    if (env_disable_mean_adjustment) {
        mean_error_adjustment_disabled = TRUE;
    }
    double base_quantile_target = 1 - base_error_quantile;

    print_config(config);
    log_info("1 - base_quantile_target: %.1e", 1-base_quantile_target);
    log_info("Disable pure base compression fallback: %s", pure_base_codec_disabled ? "TRUE" : "FALSE");



    codec_data_buffer_t codec_data_buffer;
    codec_data_buffer_init(&codec_data_buffer);

    size_t image_dims[2] = {1, config->dims[NDIMS - 1]};
    size_t tile_dims[2] = {config->dims[NDIMS - 2], config->dims[NDIMS - 1]};

    for (size_t dim_idx = 0; dim_idx < NDIMS - 1; ++dim_idx) {
        image_dims[0] *= config->dims[dim_idx];
    }

    size_t n_tiles = image_dims[0] / tile_dims[0];
    size_t tile_size = tile_dims[0] * tile_dims[1];

    // find maxval and minval
    size_t tot_size = n_tiles * tile_size;
    check_nan_inf(data, tot_size);

    float minval, maxval;
    findMinMaxf(data, tot_size, &minval, &maxval);

    int const_field = minval == maxval;

    log_trace("minval: %e, maxval: %e, const_field: %s", minval, maxval, const_field ? "TRUE" : "FALSE");

    uint16_t *scaled_data;

    if (!const_field) {
        // scale data to uint16_t
        scaled_data = (uint16_t *) malloc(tot_size * sizeof(uint16_t));
        for (size_t i = 0; i < tot_size; ++i) {
            scaled_data[i] = ((data[i] - minval) / (maxval - minval)) * (uint16_t)-1;
        }

        // encode using jpeg2000
        codec_data_buffer_clear(&codec_data_buffer);
        j2k_encode_internal(scaled_data, image_dims, tile_dims, config->base_cr, &codec_data_buffer);
        if (config->residual_compression_type == NONE) {
            base_codec_buffer = (uint8_t *) malloc(codec_data_buffer.length);
            memcpy(base_codec_buffer, codec_data_buffer.buffer, codec_data_buffer.length);
            base_codec_buffer_length = codec_data_buffer.length;
        }
        codec_data_buffer_rewind(&codec_data_buffer);
    }

    if ((config->residual_compression_type != NONE) && !const_field) {
        // decode back the image
        float *residual = (float *) malloc(tot_size * sizeof(float));
        float *residual_norm = (float *) malloc(tot_size * sizeof(float));
        float *decoded = (float *) malloc(tot_size * sizeof(float));
        j2k_decode_internal(&decoded, NULL, NULL, minval, maxval, &codec_data_buffer);

        cur_mean_error = get_mean_error(data, decoded, NULL, tot_size);

        // compute error from jpeg2000
        for (size_t i = 0; i < tot_size; ++i) {
            residual[i] = data[i] - decoded[i];
        }

        findMinMaxf(residual, tot_size, &residual_minval, &residual_maxval);
        for (size_t i = 0; i < tot_size; ++i) {
            residual_norm[i] = (residual[i] - residual_minval) / (residual_maxval - residual_minval);
        }

        if (config->residual_compression_type == MAX_ERROR ||
                config->residual_compression_type == RELATIVE_ERROR) {
            error_target = config->error, current_cr = config->base_cr;
            if (config->residual_compression_type == RELATIVE_ERROR) {
                error_target *= get_data_range(data, tot_size);
            }

            current_cr = error_bound_j2k_compression(scaled_data, image_dims, tile_dims, current_cr, &codec_data_buffer, &decoded, minval, maxval, data, tot_size, error_target, base_quantile_target);
            
            for (size_t i = 0; i < tot_size; ++i) {
                residual[i] = data[i] - decoded[i];
            }
            findMinMaxf(residual, tot_size, &residual_minval, &residual_maxval);

            float cur_max_error = fmaxf(fabsf(residual_minval), fabsf(residual_maxval));
            float best_feasible_error = -1;
            int skip_residual = cur_max_error <= error_target;
            pure_base_codec_done = base_quantile_target == 1.0;

            if (pure_base_codec_done) log_info("Pure base compression is feasible, compression error: %f, cr: %f", cur_max_error, current_cr);
            

            if (!skip_residual) {
                for (size_t i = 0; i < tot_size; ++i) {
                    residual_norm[i] = (residual[i] - residual_minval) / (residual_maxval - residual_minval);
                }
                coeffs_trunc_bits = codec_data_buffer.length * 8;
                spiht_encode(residual_norm, image_dims[0], image_dims[1], &coeffs_buf, &coeffs_size_orig, coeffs_trunc_bits, WAVELET_LEVELS);
                spiht_decode(coeffs_buf, coeffs_size_orig, residual_norm, image_dims[0], image_dims[1], coeffs_size_orig*8);
                coeffs_size = coeffs_size_orig;
                for (size_t i = 0; i < tot_size; ++i) {
                    residual[i] = residual_norm[i] * (residual_maxval - residual_minval) + residual_minval;
                }
                cur_max_error = get_max_error(data, decoded, residual, tot_size);
                if (cur_max_error > error_target) {
                    log_info("Could not reach error target of %f (%f instead), base compression max error: %f. Retry with pure base compression.", error_target, cur_max_error, fmaxf(fabsf(residual_minval), fabsf(residual_maxval)));
                    skip_residual = TRUE;
                    pure_base_codec_required = TRUE;
                    /*DONE: if this happen, go for full jpeg2000*/
                } else {
                    best_feasible_error = cur_max_error;
                    cur_mean_error = get_mean_error(data, decoded, residual, tot_size);
                }
            }
            if (!skip_residual) {
                trunc_hi = (double) coeffs_size * 8;
                trunc_lo = 112.0; // that's the number of bits in the header
                coeffs_trunc_bits = (size_t) trunc_lo;
                cur_max_error = fmaxf(fabsf(residual_minval), fabsf(residual_maxval));

                log_trace("trunc_lo: %.1f, trunc_hi: %.1f, coeffs_trunc_bytes: %lu, cur_max_error: %f, error_target: %f", trunc_lo, trunc_hi, coeffs_trunc_bits / 8, cur_max_error, error_target);

                /* TODO: scan from small values, recursive doubling*/
                
                /* TODO: exit after 64 trials or examine initial trunc_hi satisfy the error requirement*/
                best_feasible_trunc = trunc_hi;
                while (((error_target - best_feasible_error)/error_target > eps) && (trunc_hi - trunc_lo > 8 * 4)) {
                    coeffs_trunc_bits = ((size_t) ceill((trunc_hi + trunc_lo) / 2 / 8)) * 8; /* ceil to bytes*/
                    spiht_decode(coeffs_buf, coeffs_trunc_bits / 8, residual_norm, image_dims[0], image_dims[1], coeffs_trunc_bits);
                    for (size_t i = 0; i < tot_size; ++i) {
                        residual[i] = residual_norm[i] * (residual_maxval - residual_minval) + residual_minval;
                    }
                    cur_max_error = get_max_error(data, decoded, residual, tot_size);
                    if (cur_max_error > error_target) {
                        trunc_lo = coeffs_trunc_bits;
                    } else {
                        trunc_hi = coeffs_trunc_bits;
                        if (cur_max_error >= best_feasible_error) {
                            best_feasible_error = cur_max_error;
                            best_feasible_trunc = coeffs_trunc_bits;
                            cur_mean_error = get_mean_error(data, decoded, residual, tot_size);
                        }
                    }
                    log_trace("trunc_lo: %.1f, trunc_hi: %.1f, coeffs_trunc_bytes: %lu, cur_max_error: %f", trunc_lo, trunc_hi, coeffs_trunc_bits / 8, cur_max_error);
                }
                coeffs_size = (size_t) (best_feasible_trunc / 8.);
#ifdef DEBUG
                spiht_decode(coeffs_buf, coeffs_size, residual_norm, image_dims[0], image_dims[1], coeffs_size*8);
                for (size_t i = 0; i < tot_size; ++i) {
                    residual[i] = residual_norm[i] * (residual_maxval - residual_minval) + residual_minval;
                }
                cur_max_error = get_max_error(data, decoded, residual, tot_size);
                log_trace("best feasible trunc: %.1f, best feasible error: %f, actual error: %f, error_target: %f", best_feasible_trunc, best_feasible_error, cur_max_error, error_target);
                fflush(stdout);
#endif
            /* TODO: check if pure JPEG at this CR works better*/
            }
            
        }

        if (coeffs_size <= 16) coeffs_size = 0;
        
        if (coeffs_size > 0) {
            compressed_size = ZSTD_compressBound(coeffs_size * sizeof(uint8_t));
            compressed_coefficients = (uint8_t *) malloc(compressed_size);
            compressed_size = ZSTD_compress(compressed_coefficients, compressed_size, coeffs_buf, coeffs_size * sizeof(uint8_t), 22);
        }

        /* Try again with pure j2k compression , to see if adding residual compression has higher compression ratio */
        base_codec_buffer_length = codec_data_buffer.length;
        size_t base_codec_buffer_size_limit = 2 * (compressed_size + base_codec_buffer_length);
        base_codec_buffer = (uint8_t *) malloc(base_codec_buffer_size_limit);
        memcpy(base_codec_buffer, codec_data_buffer.buffer, base_codec_buffer_length);
        if ((! pure_base_codec_done) && (! pure_base_codec_disabled) && (config->residual_compression_type == MAX_ERROR ||
            config->residual_compression_type == RELATIVE_ERROR)) {
            assert(error_target > 0);
            /* ===========Maintain consistency with quantile = 0 (Not necessary) =========== */
            if (!pure_base_codec_consistency_disabled) {
                codec_data_buffer_clear(&codec_data_buffer);
                j2k_encode_internal(scaled_data, image_dims, tile_dims, config->base_cr, &codec_data_buffer);
                codec_data_buffer_rewind(&codec_data_buffer);
                j2k_decode_internal(&decoded, NULL, NULL, minval, maxval, &codec_data_buffer);
                current_cr = config->base_cr;
            }
            /* ===========Maintain consistency with quantile = 0 (Not necessary) =========== */
            error_bound_j2k_compression(scaled_data, image_dims, tile_dims, current_cr, &codec_data_buffer, &decoded, minval, maxval, data, tot_size, error_target, 1.0);

            if ((codec_data_buffer.length < compressed_size + base_codec_buffer_length) || pure_base_codec_required) {
                /* Pure JP2 is better than JP2 + SPWV */
                if (codec_data_buffer.length < compressed_size + base_codec_buffer_length)
                    log_info("Pure base compression (%lu) is better than base (%lu) + residual (%lu) compression (sum: %lu)", codec_data_buffer.length, base_codec_buffer_length, compressed_size, compressed_size + base_codec_buffer_length);
                
                cur_mean_error = get_mean_error(data, decoded, NULL, tot_size);

                compressed_size = 0;
                coeffs_size = 0;
                if (codec_data_buffer.length > base_codec_buffer_size_limit) { /* This can happen when pure_base_codec_required enabled */
                    free(base_codec_buffer);
                    base_codec_buffer = (uint8_t *) malloc(codec_data_buffer.length);
                }
                base_codec_buffer_length = codec_data_buffer.length;
                memcpy(base_codec_buffer, codec_data_buffer.buffer, base_codec_buffer_length);
            }
        }
        free(coeffs_buf);
        free(residual);
        free(residual_norm);
        free(decoded);
    }

    if (!const_field) free(scaled_data);

    log_info("Mean of compression error: %e", cur_mean_error);
    if (!mean_error_adjustment_disabled && fabs(cur_mean_error) > 1e-18) {
        minval += cur_mean_error;
        maxval += cur_mean_error;
        log_info("Adjusting minval and maxval to %f, %f", minval, maxval);
    }

    size_t codec_size = (const_field) ? sizeof(uint64_t) : base_codec_buffer_length; /*Only output array length if having a constant field*/

    size_t out_size =
            sizeof(ebcc_header_t) +
            compressed_size + codec_size;
    *out_buffer = (uint8_t *) malloc(out_size);

    log_info("coeffs_size: %lu, compressed_size: %lu, jp2_length: %lu, compression ratio: %f", coeffs_size, compressed_size, codec_size, (double) (tot_size * sizeof(float)) / out_size);

    uint8_t *iter = *out_buffer;
    ebcc_header_t header = {0};
    memcpy(header.magic, EBCC_HEADER_MAGIC, 4);
    header.version = EBCC_HEADER_VERSION;
    if (const_field) {
        header.flags |= EBCC_HEADER_FLAG_CONST_FIELD;
    }
    header.minval_bits = float_to_bits(minval);
    header.maxval_bits = float_to_bits(maxval);
    header.coeffs_size = (uint64_t) coeffs_size;
    header.residual_minval_bits = float_to_bits(residual_minval);
    header.residual_maxval_bits = float_to_bits(residual_maxval);
    header.compressed_size = (uint64_t) compressed_size;
    header.tail_size = (uint64_t) codec_size;
    memcpy(iter, &header, sizeof(header));
    iter += sizeof(header);
    if (compressed_size > 0) {
        memcpy(iter, compressed_coefficients, compressed_size);
        iter += compressed_size;
    }
    if (const_field) {
        uint64_t tot_size_u64 = (uint64_t) tot_size;
        memcpy(iter, &tot_size_u64, sizeof(uint64_t));
        iter += sizeof(uint64_t);
    } else {
        memcpy(iter, base_codec_buffer, base_codec_buffer_length);
        iter += base_codec_buffer_length;
    }
    assert((size_t) (iter - *out_buffer) == out_size);

    if (compressed_coefficients) free(compressed_coefficients);
    if (base_codec_buffer) free(base_codec_buffer);

    codec_data_buffer_destroy(&codec_data_buffer);

#ifdef ENABLE_PERF
    prctl(PR_TASK_PERF_EVENTS_DISABLE);
#endif
    return out_size;
}

size_t ebcc_encode_chunking(float *data, codec_config_t *config, uint8_t **out_buffer) {
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
    if (!dims_are_valid(chunk_dims)) {
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
        compat_config.chunk_dims[0] = 1;
        compat_config.chunk_dims[1] = compat_config.dims[1] > EBCC_MAX_INTERNAL_IMAGE_DIM ?
                1024 : compat_config.dims[1];
        compat_config.chunk_dims[2] = compat_config.dims[2] > EBCC_MAX_INTERNAL_IMAGE_DIM ?
                1024 : compat_config.dims[2];
        log_info("ebcc_encode_chunking_compat chunk dimensions: (%lu, %lu, %lu)",
                compat_config.chunk_dims[0], compat_config.chunk_dims[1], compat_config.chunk_dims[2]);
    }

    if (compat_config.residual_compression_type == RELATIVE_ERROR) {
        size_t total_size = 0;
        if (!product_size_t(compat_config.dims, &total_size) || total_size == 0) {
            log_fatal("Invalid EBCC dimensions: size overflow or zero-sized data");
            return 0;
        }
        check_nan_inf(data, total_size);
        compat_config.error *= get_data_range(data, total_size);
        compat_config.residual_compression_type = MAX_ERROR;
    }

    return ebcc_encode_chunking(data, &compat_config, out_buffer);
}

void j2k_decode_internal(float **data, size_t *height, size_t *width, float minval, float maxval,
        codec_data_buffer_t *codec_data_buffer) {
    codec_data_buffer_rewind(codec_data_buffer);

    opj_stream_t *stream = opj_stream_default_create(OPJ_TRUE);

    opj_stream_set_user_data(stream, codec_data_buffer, NULL);
    opj_stream_set_user_data_length(stream, codec_data_buffer->length);
    opj_stream_set_read_function(stream, (opj_stream_read_fn) read_from_buffer_stream);

    opj_dparameters_t param;
    opj_set_default_decoder_parameters(&param);

    param.decod_format = 0;
    param.cp_layer = 0;
    param.cp_reduce = 0;

    opj_codec_t *codec = opj_create_decompress(OPJ_CODEC_J2K);

    opj_setup_decoder(codec, &param);

    opj_image_t *image;
    opj_read_header(stream, codec, &image);

    opj_decode(codec, stream, image);
    opj_end_decompress(codec, stream);

    if (width) {
        *width = image->x1 - image->x0;
    }
    if (height) {
        *height = image->y1 - image->y0;
    }
    size_t num_pixels = (image->x1 - image->x0) * (image->y1 - image->y0);
    if (!*data) {
        *data = (float *) malloc(num_pixels * sizeof(float));
    }
    for (size_t i = 0; i < num_pixels; ++i) {
        (*data)[i] = ((float) image->comps[0].data[i] / (uint16_t)-1) * (maxval - minval) + minval;
    }

    opj_stream_destroy(stream);
    opj_destroy_codec(codec);
    opj_image_destroy(image);
}

static int legacy_read_bytes(const uint8_t **iter, const uint8_t *end, void *dst, size_t nbytes) {
    if ((size_t) (end - *iter) < nbytes) {
        return FALSE;
    }
    memcpy(dst, *iter, nbytes);
    *iter += nbytes;
    return TRUE;
}

size_t ebcc_decode_legacy(uint8_t *data, size_t data_size, float **out_buffer) {
    codec_data_buffer_t codec_data_buffer;
    size_t tot_size = 0;
    const uint8_t *iter = data;
    const uint8_t *end = data + data_size;

    float minval = 0.0f, maxval = 0.0f, residual_minval = 0.0f, residual_maxval = 0.0f;
    size_t coeffs_size = 0, compressed_coefficient_size = 0;
    if (!legacy_read_bytes(&iter, end, &minval, sizeof(float)) ||
            !legacy_read_bytes(&iter, end, &maxval, sizeof(float)) ||
            !legacy_read_bytes(&iter, end, &coeffs_size, sizeof(size_t)) ||
            !legacy_read_bytes(&iter, end, &residual_minval, sizeof(float)) ||
            !legacy_read_bytes(&iter, end, &residual_maxval, sizeof(float)) ||
            !legacy_read_bytes(&iter, end, &compressed_coefficient_size, sizeof(size_t))) {
        log_fatal("Invalid legacy encoded data: truncated header");
        return 0;
    }

    if ((size_t) (end - iter) < compressed_coefficient_size) {
        log_fatal("Invalid legacy encoded data: truncated residual payload");
        return 0;
    }
    uint8_t *coefficient_data = (uint8_t *) iter;
    iter += compressed_coefficient_size;

    size_t height = 0, width = 0;
    int const_field = minval == maxval;
    if (const_field) {
        if (!legacy_read_bytes(&iter, end, &tot_size, sizeof(size_t))) {
            log_fatal("Invalid legacy encoded data: missing const-field length");
            return 0;
        }
        *out_buffer = (float *) malloc(tot_size * sizeof(float));
        for (size_t i = 0; i < tot_size; ++i) {
            (*out_buffer)[i] = minval;
        }
    } else {
        codec_data_buffer.buffer = (uint8_t *) iter;
        codec_data_buffer.size = (size_t) (end - iter);
        codec_data_buffer.length = (size_t) (end - iter);
        codec_data_buffer.offset = 0;
        j2k_decode_internal(out_buffer, &height, &width, minval, maxval, &codec_data_buffer);
        tot_size = height * width;
    }

    log_info("Compression Ratio: %f", (double) (tot_size * sizeof(float)) / data_size);

    if (compressed_coefficient_size > 0 && coeffs_size > 0) {
        if (const_field) {
            log_fatal("Invalid legacy encoded data: residual data cannot be applied to const field");
            return 0;
        }
        uint8_t *coeffs = (uint8_t *) calloc(coeffs_size, sizeof(uint8_t));
        float *residual = (float *) calloc(tot_size, sizeof(float));
        ZSTD_decompress(coeffs, coeffs_size * sizeof(uint8_t), coefficient_data, compressed_coefficient_size);
        spiht_decode(coeffs, coeffs_size, residual, height, width, coeffs_size * 8);

        for (size_t i = 0; i < tot_size; ++i) {
            (*out_buffer)[i] += residual[i] * (residual_maxval - residual_minval) + residual_minval;
        }

        free(coeffs);
        free(residual);
    }

    return tot_size;
}

size_t ebcc_decode(uint8_t *data, size_t data_size, float **out_buffer) {
    assert_endianness();

    codec_data_buffer_t codec_data_buffer;
    size_t tot_size = 0;
    uint8_t *iter = data;

    if (data_size < sizeof(ebcc_header_t) || memcmp(data, EBCC_HEADER_MAGIC, 4) != 0) {
        return ebcc_decode_legacy(data, data_size, out_buffer);
    }

    ebcc_header_t header;
    memcpy(&header, iter, sizeof(header));
    iter += sizeof(header);

    if (header.version != EBCC_HEADER_VERSION) {
        log_fatal("Unsupported EBCC header version: %u", header.version);
        return 0;
    }

    if (header.coeffs_size > SIZE_MAX || header.compressed_size > SIZE_MAX || header.tail_size > SIZE_MAX) {
        log_fatal("Invalid encoded data: header field exceeds local SIZE_MAX");
        return 0;
    }

    float minval = bits_to_float(header.minval_bits);
    float maxval = bits_to_float(header.maxval_bits);
    float residual_minval = bits_to_float(header.residual_minval_bits);
    float residual_maxval = bits_to_float(header.residual_maxval_bits);
    size_t coeffs_size = (size_t) header.coeffs_size;
    size_t compressed_coefficient_size = (size_t) header.compressed_size;
    size_t tail_size = (size_t) header.tail_size;

    size_t header_and_payload_size = sizeof(ebcc_header_t);
    if (compressed_coefficient_size > data_size - header_and_payload_size) {
        log_fatal("Invalid encoded data: truncated payload");
        return 0;
    }
    header_and_payload_size += compressed_coefficient_size;
    if (tail_size > data_size - header_and_payload_size) {
        log_fatal("Invalid encoded data: truncated payload");
        return 0;
    }
    header_and_payload_size += tail_size;

    uint8_t *coefficient_data = iter;
    iter += compressed_coefficient_size;

    size_t height = 0, width = 0;
    int const_field = (header.flags & EBCC_HEADER_FLAG_CONST_FIELD) != 0;
    if (const_field) {
        if (tail_size != sizeof(uint64_t)) {
            log_fatal("Invalid encoded data: const-field payload must contain uint64_t length");
            return 0;
        }
        uint64_t tot_size_u64 = 0;
        memcpy(&tot_size_u64, iter, sizeof(uint64_t));
        if (tot_size_u64 > SIZE_MAX) {
            log_fatal("Invalid encoded data: const-field size exceeds local SIZE_MAX");
            return 0;
        }
        tot_size = (size_t) tot_size_u64;
        iter += sizeof(uint64_t);
        *out_buffer = (float *) malloc(tot_size * sizeof(float));
        for (size_t i = 0; i < tot_size; ++i) {
            (*out_buffer)[i] = minval;
        }
    } else {
        codec_data_buffer.buffer = iter;
        codec_data_buffer.size = tail_size;
        codec_data_buffer.length = tail_size;
        codec_data_buffer.offset = 0;
        j2k_decode_internal(out_buffer, &height, &width, minval, maxval, &codec_data_buffer);
        tot_size = height * width;
        iter += tail_size;
    }

    log_info("Compression Ratio: %f", (double) (tot_size * sizeof(float)) / data_size);

    if (compressed_coefficient_size > 0 && coeffs_size > 0) {
        if (const_field) {
            log_fatal("Invalid encoded data: residual data cannot be applied to const field");
            return 0;
        }
        uint8_t *coeffs = (uint8_t *) calloc(coeffs_size, sizeof(uint8_t));
        float *residual = (float *) calloc(tot_size, sizeof(float));
        ZSTD_decompress(coeffs, coeffs_size * sizeof(uint8_t), coefficient_data, compressed_coefficient_size);


        spiht_decode(coeffs, coeffs_size, residual, height, width, coeffs_size * 8);

        for (size_t i = 0; i < tot_size; ++i) {
            (*out_buffer)[i] += residual[i] * (residual_maxval - residual_minval) + residual_minval;
        }

        free(coeffs);
        free(residual);
    }

    if ((size_t) (iter - data) != header_and_payload_size) {
        log_fatal("Invalid encoded data: payload size mismatch");
        return 0;
    }

    return tot_size;
}

size_t ebcc_decode_chunking(uint8_t *data, size_t data_size, float **out_buffer) {
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
    if (!dims_are_valid(chunk_dims)) {
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

    *out_buffer = (float *) malloc(total_bytes);
    if (!*out_buffer) {
        log_fatal("Failed to allocate chunked EBCC decode output");
        return 0;
    }

    for (size_t chunk_linear = 0; chunk_linear < header_num_chunks; ++chunk_linear) {
        if ((size_t) (end - iter) < sizeof(uint64_t)) {
            log_fatal("Invalid chunked EBCC data: missing chunk size");
            free(*out_buffer);
            *out_buffer = NULL;
            return 0;
        }

        uint64_t chunk_stream_size_u64 = 0;
        memcpy(&chunk_stream_size_u64, iter, sizeof(chunk_stream_size_u64));
        iter += sizeof(chunk_stream_size_u64);
        size_t chunk_stream_size = 0;
        if (!uint64_to_size_t(chunk_stream_size_u64, &chunk_stream_size)) {
            log_fatal("Invalid chunked EBCC data: chunk payload size exceeds local SIZE_MAX");
            free(*out_buffer);
            *out_buffer = NULL;
            return 0;
        }
        if (chunk_stream_size > (size_t) (end - iter)) {
            log_fatal("Invalid chunked EBCC data: truncated chunk payload");
            free(*out_buffer);
            *out_buffer = NULL;
            return 0;
        }

        float *chunk_buffer = NULL;
        size_t decoded_chunk_size = ebcc_decode((uint8_t *) iter, chunk_stream_size, &chunk_buffer);
        if (decoded_chunk_size != header_chunk_size || !chunk_buffer) {
            log_fatal("Invalid chunked EBCC data: decoded chunk %lu has %lu values, expected %lu",
                    chunk_linear, decoded_chunk_size, header_chunk_size);
            free_buffer(chunk_buffer);
            free(*out_buffer);
            *out_buffer = NULL;
            return 0;
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
        free(*out_buffer);
        *out_buffer = NULL;
        return 0;
    }

    return total_size;
}

void free_buffer(void *buffer) {
    if (buffer) {
        free(buffer);
    }
}
