#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "log.h"

#include "openjpeg.h"
#include "zstd.h"
#include <assert.h>
#include <malloc.h>
#include "spiht_re.h"


#define MIN(x, y) ((x)<(y)?(x):(y))
#define CEIL(x,y) (((x) + (y) - 1) / (y))
#define TRUE 1
#define FALSE 0

#define WAVELET_LEVELS 3

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

void codec_data_buffer_reset(codec_data_buffer_t* data) {
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

typedef enum {
    NONE,
    SPARSIFICATION_FACTOR,
    MAX_ERROR,
    RELATIVE_ERROR,
    QUANTILE
} residual_t;

const char* residual_t_names[] = {
    "NONE",
    "SPARSIFICATION_FACTOR",
    "MAX_ERROR",
    "RELATIVE_ERROR",
    "QUANTILE"
};

#define NDIMS 3

typedef struct {
    size_t dims[NDIMS];
    float base_cr;
    residual_t residual_compression_type;
    float residual_cr;
    float error;
    double quantile;
} codec_config_t;

void print_config(codec_config_t *config) {
    log_info("dimensions:\t(%lu, %lu, %lu)", config->dims[0], config->dims[1], config->dims[2]);
    log_info("base_cr:\t%f", config->base_cr);
    log_info("residual type:\t%s", residual_t_names[config->residual_compression_type]);
    switch (config->residual_compression_type) {
        case NONE:
            break;
        case SPARSIFICATION_FACTOR:
            log_info("sparsification:\t%f", config->residual_cr);
            break;
        case MAX_ERROR:
            log_info("max error:\t%f", config->error);
            break;
        case RELATIVE_ERROR:
            log_info("relative error:\t%f", config->error);
            break;
        case QUANTILE:
            log_info("quantile:\t%f", config->quantile);
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
    codec_data_buffer_init(codec_data_buffer);
    j2k_encode_internal(scaled_data, image_dims, tile_dims, current_cr, codec_data_buffer);
    codec_data_buffer_reset(codec_data_buffer);
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

size_t encode_climate_variable(float *data, codec_config_t *config, uint8_t **out_buffer) {
    int pure_j2k_required = FALSE, pure_j2k_done = FALSE, pure_j2k_disabled = FALSE, pure_j2k_consistency_disabled = FALSE;
    size_t compressed_size = 0, jp2_buffer_length = 0;
    uint8_t *compressed_coefficients = NULL;
    uint8_t *coeffs_buf = NULL;
    uint8_t *jp2_buffer = NULL; 
    float residual_maxval = 0., residual_minval = 0., error_target = -1, current_cr = -1;
    size_t coeffs_size = 0, coeffs_size_orig = 0, coeffs_trunc_bits = 0; /*coeffs_size: #bytes*/
    double trunc_hi, trunc_lo, best_feasible_trunc;
    double quantile = config->quantile, eps = 1e-8, base_error_quantile=1e-6;

    // Load base_error_quantile from env var EBCC_INIT_BASE_ERROR_QUANTILE, default to 1e-6 if not set

    log_set_level_from_env();
    
    const char *env_base_error_quantile = getenv("EBCC_INIT_BASE_ERROR_QUANTILE");
    const char *env_disable_pure_jp2_fallback = getenv("EBCC_DISABLE_PURE_JP2_FALLBACK");
    const char *env_disable_pure_jp2_consistency = getenv("EBCC_DISABLE_PURE_JP2_FALLBACK_CONSISTENCY");
    if (env_base_error_quantile) {
        base_error_quantile = strtod(env_base_error_quantile, NULL);
    }
    if (env_disable_pure_jp2_fallback) {
        pure_j2k_disabled = TRUE;
    }
    if (env_disable_pure_jp2_consistency) {
        pure_j2k_consistency_disabled = TRUE;
    }
    double base_quantile_target = 1 - base_error_quantile;

    print_config(config);
    log_info("1 - base_quantile_target: %.1e", 1-base_quantile_target);
    log_info("Disable pure JP2 fallback: %s", pure_j2k_disabled ? "TRUE" : "FALSE");



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
        j2k_encode_internal(scaled_data, image_dims, tile_dims, config->base_cr, &codec_data_buffer);
        if (config->residual_compression_type == NONE) {
            memcpy(jp2_buffer, codec_data_buffer.buffer, codec_data_buffer.length);
            jp2_buffer_length = codec_data_buffer.length;
        }
        codec_data_buffer_reset(&codec_data_buffer);
    }

    if ((config->residual_compression_type != NONE) && !const_field) {
        // decode back the image
        float *residual = (float *) malloc(tot_size * sizeof(float));
        float *residual_norm = (float *) malloc(tot_size * sizeof(float));
        float *decoded = (float *) malloc(tot_size * sizeof(float));
        j2k_decode_internal(&decoded, NULL, NULL, minval, maxval, &codec_data_buffer);

        // compute error from jpeg2000
        for (size_t i = 0; i < tot_size; ++i) {
            residual[i] = data[i] - decoded[i];
        }

        findMinMaxf(residual, tot_size, &residual_minval, &residual_maxval);
        for (size_t i = 0; i < tot_size; ++i) {
            residual_norm[i] = (residual[i] - residual_minval) / (residual_maxval - residual_minval);
        }

        if (config->residual_compression_type == QUANTILE) {
            assert(0); /* Deprecated */
        } else if (config->residual_compression_type == SPARSIFICATION_FACTOR) {
            
            double q_ratio = 1. - (1. / config->residual_cr);
            coeffs_trunc_bits = tot_size * 8 - (size_t) ((double)tot_size / config->residual_cr * 8 );
            spiht_encode(residual_norm, image_dims[0], image_dims[1], &coeffs_buf, &coeffs_size_orig, coeffs_trunc_bits, WAVELET_LEVELS);
            coeffs_size = coeffs_size_orig;
        } else if (config->residual_compression_type == MAX_ERROR ||
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
            pure_j2k_done = base_quantile_target == 1.0;

            if (pure_j2k_done) log_info("Pure JP2 compression is feasible, compression error: %f, cr: %f", cur_max_error, current_cr);
            

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
                    log_warn("Could not reach error target of %f (%f instead), base compression max error: %f. Retry with pure base compression.", error_target, cur_max_error, fmaxf(fabsf(residual_minval), fabsf(residual_maxval)));
                    skip_residual = TRUE;
                    pure_j2k_required = TRUE;
                    /*DONE: if this happen, go for full jpeg2000*/
                } else {
                    best_feasible_error = cur_max_error;
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

        log_info("coeffs_size: %lu, compressed_size: %lu, jp2_length: %lu, compression ratio: %f", coeffs_size, compressed_size, codec_data_buffer.length, (float) (tot_size * sizeof(float)) / (compressed_size + codec_data_buffer.length));

        /* Try again with pure j2k compression , to see if adding residual compression has higher compression ratio */
        jp2_buffer_length = codec_data_buffer.length;
        size_t jp2_buffer_size_limit = 2 * (compressed_size + jp2_buffer_length);
        jp2_buffer = (uint8_t *) malloc(jp2_buffer_size_limit);
        memcpy(jp2_buffer, codec_data_buffer.buffer, jp2_buffer_length);
        if ((! pure_j2k_done) && (! pure_j2k_disabled) && (config->residual_compression_type == MAX_ERROR ||
            config->residual_compression_type == RELATIVE_ERROR)) {
            assert(error_target > 0);
            /* ===========Maintain consistency with quantile = 0 (Not necessary) =========== */
            if (!pure_j2k_consistency_disabled) {
                codec_data_buffer_init(&codec_data_buffer);
                j2k_encode_internal(scaled_data, image_dims, tile_dims, config->base_cr, &codec_data_buffer);
                codec_data_buffer_reset(&codec_data_buffer);
                j2k_decode_internal(&decoded, NULL, NULL, minval, maxval, &codec_data_buffer);
                current_cr = config->base_cr;
            }
            /* ===========Maintain consistency with quantile = 0 (Not necessary) =========== */
            error_bound_j2k_compression(scaled_data, image_dims, tile_dims, current_cr, &codec_data_buffer, &decoded, minval, maxval, data, tot_size, error_target, 1.0);

            if ((codec_data_buffer.length < compressed_size + jp2_buffer_length) || pure_j2k_required) {
                /* Pure JP2 is better than JP2 + SPWV */
                log_info("Pure JP2 (%lu) is better than JP2 (%lu) + SPWV (%lu) (sum: %lu)", codec_data_buffer.length, jp2_buffer_length, compressed_size, compressed_size + jp2_buffer_length);

                compressed_size = 0;
                coeffs_size = 0;
                if (codec_data_buffer.length > jp2_buffer_size_limit) { /* This can happen when pure_j2k_required enabled */
                    free(jp2_buffer);
                    jp2_buffer = (uint8_t *) malloc(codec_data_buffer.length);
                }
                jp2_buffer_length = codec_data_buffer.length;
                memcpy(jp2_buffer, codec_data_buffer.buffer, jp2_buffer_length);
            }
        }
        free(coeffs_buf);
        free(residual);
        free(residual_norm);
        free(decoded);
    }

    if (!const_field) free(scaled_data);

    size_t codec_size = (const_field) ? sizeof(size_t) : jp2_buffer_length; /*Only output array length if having a constant field*/

    size_t out_size =
            2 * sizeof(float) /* minval and maxval */ +
            sizeof(size_t) + /* coeffs size */
            2 * sizeof(float) /* residual_minval, residual_maxval */ +
            sizeof(size_t) /* compressed_size */ + 
            compressed_size + codec_size;
    *out_buffer = (uint8_t *) malloc(out_size);

    uint8_t *iter = *out_buffer;
    memcpy(iter, &minval, sizeof(float));
    iter += sizeof(float);
    memcpy(iter, &maxval, sizeof(float));
    iter += sizeof(float);
    memcpy(iter, &coeffs_size, sizeof(size_t));
    iter += sizeof(size_t);
    memcpy(iter, &residual_minval, sizeof(float));
    iter += sizeof(float);
    memcpy(iter, &residual_maxval, sizeof(float));
    iter += sizeof(float);
    memcpy(iter, &compressed_size, sizeof(size_t));
    iter += sizeof(size_t);
    if (compressed_size > 0) {
        memcpy(iter, compressed_coefficients, compressed_size);
        iter += compressed_size;
    }
    if (const_field) {
        memcpy(iter, &tot_size, sizeof(size_t));
    } else {
        memcpy(iter, jp2_buffer, jp2_buffer_length);
    }
    assert(iter - *out_buffer == out_size - codec_size);

    if (compressed_coefficients) free(compressed_coefficients);
    free(jp2_buffer);

    codec_data_buffer_destroy(&codec_data_buffer);

    return out_size;
}

void j2k_decode_internal(float **data, size_t *height, size_t *width, float minval, float maxval,
        codec_data_buffer_t *codec_data_buffer) {
    codec_data_buffer_reset(codec_data_buffer);

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

size_t decode_climate_variable(uint8_t *data, size_t data_size, float **out_buffer) {
    codec_data_buffer_t codec_data_buffer;
    size_t tot_size = 0;
    uint8_t *iter = data;
    float minval = *((float *) iter);
    iter += sizeof(float);
    float maxval = *((float *) iter);
    iter += sizeof(float);
    size_t coeffs_size = *((size_t *) iter);
    iter += sizeof(size_t);
    float residual_minval = *((float *) iter);
    iter += sizeof(float);
    float residual_maxval = *((float *) iter);
    iter += sizeof(float);
    size_t compressed_coefficient_size = *((size_t *) iter);
    iter += sizeof(size_t);
    uint8_t *coefficient_data = iter;
    iter += compressed_coefficient_size;

    size_t height, width;
    int const_field = minval == maxval;
    if (const_field) {
        tot_size = *((size_t *) iter);
        iter += sizeof(size_t);
        *out_buffer = (float *) malloc(tot_size * sizeof(float));
        for (size_t i = 0; i < tot_size; ++i) {
            (*out_buffer)[i] = minval;
        }
    } else {
        codec_data_buffer.buffer = iter;
        codec_data_buffer.size = data_size - (iter - data);
        codec_data_buffer.length = data_size - (iter - data);
        codec_data_buffer.offset = 0;
        j2k_decode_internal(out_buffer, &height, &width, minval, maxval, &codec_data_buffer);
        tot_size = height * width;
    }

    if (compressed_coefficient_size > 0 && coeffs_size > 0) {
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
