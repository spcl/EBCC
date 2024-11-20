#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "openjpeg.h"
#include "wavelet.h"
#include "quantile.h"
#include "zstd.h"
#include <assert.h>
#include <malloc.h>
#include "spiht.h"

#define MIN(x, y) ((x)<(y)?(x):(y))
#define CEIL(x,y) (((x) + (y) - 1) / (y))

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
    printf("dimensions:\t(%lu, %lu, %lu)\n", config->dims[0], config->dims[1], config->dims[2]);
    printf("base_cr:\t%f\n", config->base_cr);
    printf("residual type:\t%s\n", residual_t_names[config->residual_compression_type]);
    switch (config->residual_compression_type) {
        case NONE:
            break;
        case SPARSIFICATION_FACTOR:
            printf("sparsification:\t%f\n", config->residual_cr);
            break;
        case MAX_ERROR:
            printf("max error:\t%f\n", config->error);
            break;
        case RELATIVE_ERROR:
            printf("relative error:\t%f\n", config->error);
            break;
        case QUANTILE:
            printf("quantile:\t%f\n", config->quantile);
            break;
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
        float cur_error = fabsf(data[i] - decoded[i] - residual_value);
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
        float cur_error = fabsf(data[i] - decoded[i] - residual_value);
        if (cur_error > error_target) {
            n++;
        }
    }
    return 1. -  ((double) n / tot_size);
}

#define STOP_CR 50

void findMinMaxf(const float *array, size_t size, float *min, float *max) {
    float min_val = INFINITY, max_val = -INFINITY;
    if (size == 0) {
        return;
    }

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

void sparsify_coefficients(const double *coeffs, double *coeffs_copy, const size_t coeffs_size, float *residual, const size_t image_dims[2], 
                           const float *data, const float *decoded, const codec_config_t *config, const size_t tot_size) {
    float residual_cr = 2000; /*50*/
    float cur_max_error = 0, best_max_error = 0;
    double quantile, data_range = 1.0, *coeffs_best;

    coeffs_best = (double *) malloc(coeffs_size * sizeof(double));

    memset(coeffs_copy, 0, coeffs_size * sizeof(double));
    memset(coeffs_best, 0, coeffs_size * sizeof(double));
    if (config->residual_compression_type == RELATIVE_ERROR) {
        data_range = get_data_range(data, tot_size);
        cur_max_error = get_max_relative_error(data, decoded, NULL, tot_size, data_range);
    } else {
        cur_max_error = get_max_error(data, decoded, NULL, tot_size);
    }
    best_max_error = cur_max_error;

    while (cur_max_error > config->error && residual_cr >= STOP_CR) {
        memcpy(coeffs_copy, coeffs, coeffs_size * sizeof(double));
        double q_ratio = 1. - (1. / residual_cr);
        quantile = zero_out_quantile(coeffs_copy, coeffs_size, q_ratio);
        wavelib_backward(residual, image_dims[0], image_dims[1], WAVELET_LEVELS, coeffs_copy);

        if (config->residual_compression_type == RELATIVE_ERROR) {
            cur_max_error = get_max_relative_error(data, decoded, residual, tot_size, data_range);
        } else {
            cur_max_error = get_max_error(data, decoded, residual, tot_size);
        }
        if (cur_max_error < best_max_error) {
            best_max_error = cur_max_error;
            memcpy(coeffs_best, coeffs_copy, coeffs_size * sizeof(double));
        }
#ifdef DEBUG
        printf("Current max error: %f (ABS %f), residual_cr: %f\n", cur_max_error, cur_max_error*data_range, residual_cr);
#endif
        residual_cr /= sqrtf(2.f);
    } 

    if (cur_max_error > config->error) {
        fprintf(stderr, "Could not reach error target of %f (%f instead).\n", config->error, best_max_error);
        memcpy(coeffs_copy, coeffs_best, coeffs_size * sizeof(double));
    }
    free(coeffs_best);
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
    /* TODO: log down best feasible cr!, best feasible error*/
    /* TODO: take error target quantile from env*/
    /* TODO: log according to env*/
    while (error_target_quantile < base_quantile_target) {
        cr_lo /= 2;
        error_target_quantile = emulate_j2k_compression(scaled_data, image_dims, tile_dims, cr_lo, codec_data_buffer, decoded, minval, maxval, data, tot_size, error_target);
#ifdef DEBUG
        printf("cr_lo: %f, 1-error_target_quantile: %.1e, jp2_length: %lu\n", cr_lo, 1-error_target_quantile, codec_data_buffer->length);
        fflush(stdout);
#endif
    }
    error_target_quantile = error_target_quantile_prev;
    while (error_target_quantile >= base_quantile_target) {
        cr_hi *= 2;
        error_target_quantile = emulate_j2k_compression(scaled_data, image_dims, tile_dims, cr_hi, codec_data_buffer, decoded, minval, maxval, data, tot_size, error_target);
#ifdef DEBUG
        printf("cr_hi: %f, 1-error_target_quantile: %.1e, jp2_length: %lu\n", cr_hi, 1-error_target_quantile, codec_data_buffer->length);
        fflush(stdout);
#endif
    }
    error_target_quantile = error_target_quantile_prev;

    assert(cr_lo <= cr_hi);
    while ((fabs(error_target_quantile - base_quantile_target) > eps || error_target_quantile == 1.0) && cr_hi - cr_lo > 1.) {
        current_cr = (cr_lo + cr_hi) / 2;
        error_target_quantile = emulate_j2k_compression(scaled_data, image_dims, tile_dims, current_cr, codec_data_buffer, decoded, minval, maxval, data, tot_size, error_target);
#ifdef DEBUG
        printf("current_cr: %f, 1-error_target_quantile: %.1e, jp2_length: %lu\n", current_cr, 1-error_target_quantile, codec_data_buffer->length);
        fflush(stdout);
#endif
        if (error_target_quantile < base_quantile_target) {
            cr_hi = current_cr;
        } else {
            cr_lo = current_cr;
        }
    }
    /* use cr_lo as the best feasible cr */
    error_target_quantile = emulate_j2k_compression(scaled_data, image_dims, tile_dims, cr_lo, codec_data_buffer, decoded, minval, maxval, data, tot_size, error_target);
    if (error_target_quantile < base_quantile_target) {
        fprintf(stderr, "Could not reach error target quantile of (1-%.2e) (1-%.2e instead).\n", 1-base_quantile_target, 1-error_target_quantile);
    }
    return cr_lo;

}

size_t encode_climate_variable(float *data, codec_config_t *config, uint8_t **out_buffer) {
#ifdef DEBUG
    print_config(config);
#endif
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
    float minval = data[0], maxval = data[0];
    for (size_t i = 1; i < tot_size; ++i) {
        if (data[i] < minval) {
            minval = data[i];
        } else if (data[i] > maxval) {
            maxval = data[i];
        }
    }

    // scale data to uint16_t
    uint16_t *scaled_data = (uint16_t *) malloc(tot_size * sizeof(uint16_t));
    for (size_t i = 0; i < tot_size; ++i) {
        scaled_data[i] = ((data[i] - minval) / (maxval - minval)) * (uint16_t)-1;
    }

    // encode using jpeg2000
    j2k_encode_internal(scaled_data, image_dims, tile_dims, config->base_cr, &codec_data_buffer);

    
    codec_data_buffer_reset(&codec_data_buffer);

    size_t compressed_size = 0;
    uint8_t *compressed_coefficients = NULL;
    uint8_t *coeffs_buf = NULL;
    float residual_maxval = 0., residual_minval = 0.;
    size_t coeffs_size = 0, coeffs_size_orig = 0, coeffs_trunc_bits = 0; /*coeffs_size: #bytes*/
    double trunc_hi, trunc_lo, best_feasible_trunc;
    double quantile = config->quantile, base_quantile_target = 1-1e-6, eps=1e-8;
#ifdef DEBUG
    printf("1 - base_quantile_target: %.1e\n", 1-base_quantile_target);
#endif
    if (config->residual_compression_type != NONE) {
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
            float error_target = config->error, current_cr = config->base_cr;
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
                    fprintf(stderr, "Could not reach error target of %f (%f instead), base compression max error: %f.\n Retry with pure base compression.\n", error_target, cur_max_error, fmaxf(fabsf(residual_minval), fabsf(residual_maxval)));
                    skip_residual = 1;
                    /*DONE: if this happen, go for full jpeg2000*/
                    current_cr = error_bound_j2k_compression(scaled_data, image_dims, tile_dims, current_cr*0.8, &codec_data_buffer, &decoded, minval, maxval, data, tot_size, error_target, 1.0);
                } else {
                    best_feasible_error = cur_max_error;
                }
            }
            if (!skip_residual) {
                trunc_hi = (double) coeffs_size * 8;
                trunc_lo = 112.0; // that's the number of bits in the header
                coeffs_trunc_bits = (size_t) trunc_lo;
                cur_max_error = fmaxf(fabsf(residual_minval), fabsf(residual_maxval));
#ifdef DEBUG
                printf("trunc_lo: %.1f, trunc_hi: %.1f, coeffs_trunc_bytes: %lu, cur_max_error: %f, error_target: %f\n", trunc_lo, trunc_hi, coeffs_trunc_bits / 8, cur_max_error, error_target);
                fflush(stdout);
#endif
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
#ifdef DEBUG
                    printf("trunc_lo: %.1f, trunc_hi: %.1f, coeffs_trunc_bytes: %lu, cur_max_error: %f\n", trunc_lo, trunc_hi, coeffs_trunc_bits / 8, cur_max_error);
                    fflush(stdout);
#endif
                }
                coeffs_size = (size_t) (best_feasible_trunc / 8.);
#ifdef DEBUG
                spiht_decode(coeffs_buf, coeffs_size, residual_norm, image_dims[0], image_dims[1], coeffs_size*8);
                for (size_t i = 0; i < tot_size; ++i) {
                    residual[i] = residual_norm[i] * (residual_maxval - residual_minval) + residual_minval;
                }
                cur_max_error = get_max_error(data, decoded, residual, tot_size);
                printf("best feasible trunc: %.1f, best feasible error: %f, actual error: %f, error_target: %f\n", best_feasible_trunc, best_feasible_error, cur_max_error, error_target);
                fflush(stdout);
#endif
            /* TODO: check if pure JPEG at this CR works better*/
            }
            
        }

        free(scaled_data);

        if (coeffs_size <= 16) coeffs_size = 0;
        
        if (coeffs_size > 0) {
            compressed_size = ZSTD_compressBound(coeffs_size * sizeof(uint8_t));
            compressed_coefficients = (uint8_t *) malloc(compressed_size);
            compressed_size = ZSTD_compress(compressed_coefficients, compressed_size, coeffs_buf, coeffs_size * sizeof(uint8_t), 22);
        }
#ifdef DEBUG
        printf("coeffs_size: %lu, compressed_size: %lu, jp2_length: %lu, compression ratio: %f\n", coeffs_size, compressed_size, codec_data_buffer.length, (float) (tot_size * sizeof(float)) / (compressed_size + codec_data_buffer.length));
#endif

        free(coeffs_buf);
        free(residual);
        free(residual_norm);
        free(decoded);
    }

    size_t out_size =
            2 * sizeof(float) /* minval and maxval */ +
            sizeof(size_t) + /* coeffs size */
            2 * sizeof(float) /* residual_minval, residual_maxval */ +
            sizeof(size_t) /* compressed_size */ + 
            compressed_size + codec_data_buffer.length;
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
    memcpy(iter, compressed_coefficients, compressed_size);
    iter += compressed_size;
    memcpy(iter, codec_data_buffer.buffer, codec_data_buffer.length);

    assert(iter - *out_buffer == out_size - codec_data_buffer.length);

    if (compressed_coefficients) free(compressed_coefficients);

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

    codec_data_buffer.buffer = iter;
    codec_data_buffer.size = data_size - (iter - data);
    codec_data_buffer.length = data_size - (iter - data);
    codec_data_buffer.offset = 0;

    size_t height, width;
    j2k_decode_internal(out_buffer, &height, &width, minval, maxval, &codec_data_buffer);

    size_t tot_size = height * width;

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