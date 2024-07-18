#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <openjpeg.h>

#include "wavelet.h"
#include "quantile.h"
#include "zstd.h"
#include <assert.h>
#include <malloc.h>

#define MIN(x, y) ((x)<(y)?(x):(y))
#define CEIL(x,y) (((x) + (y) - 1) / (y))

#define WAVELET_LEVELS 3

typedef struct {
    uint8_t *buffer;
    size_t size;        // size of the buffer (maximum allowed storage)
    size_t length;      // number of data bytes in the buffer
    size_t offset;
} codec_data_buffer_t;

void j2k_decode_internal(float **data, size_t *height, size_t *width, float minval, float maxval, codec_data_buffer_t *codec_data_buffer);

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

void j2k_encode_internal(void *data, size_t *image_dims, size_t *tile_dims, float base_cr, codec_data_buffer_t *codec_data_buffer) {
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
            opj_write_tile(codec, i, ((OPJ_BYTE *) data) + i * (tile_size * sizeof(uint16_t)), tile_size * sizeof(uint16_t), stream);
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
    QUANTILE
} residual_t;

const char* residual_t_names[] = {
    "NONE",
    "SPARSIFICATION_FACTOR",
    "MAX_ERROR",
    "QUANTILE"
};

#define NDIMS 3

typedef struct {
    size_t dims[NDIMS];
    float base_cr;
    residual_t residual_compression_type;
    float residual_cr;
    float max_error;
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
            printf("max error:\t%f\n", config->max_error);
            break;
        case QUANTILE:
            printf("quantile:\t%f\n", config->quantile);
            break;
    }
}

void spiht_encode(float *buffer, size_t height, size_t width, float bit_rate, void **out_buffer, size_t *out_size);
void spiht_decode(void *buffer, size_t size, float *out_buffer);

size_t encode_climate_variable(float *data, codec_config_t *config, uint8_t **out_buffer) {
    print_config(config);

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

    free(scaled_data);
    codec_data_buffer_reset(&codec_data_buffer);

    size_t compressed_size = 0;
    uint8_t *compressed_coefficients = NULL;

    double coeffs_min = 0, coeffs_max = 0;
    double *coeffs = NULL;
    size_t coeffs_size = 0;
    size_t coo_size = 0;
    double quantile = config->quantile;
    if (config->residual_compression_type != NONE) {
        // decode back the image
        float *residual = (float *) malloc(tot_size * sizeof(float));
        float *decoded = (float *) malloc(tot_size * sizeof(float));
        j2k_decode_internal(&decoded, NULL, NULL, minval, maxval, &codec_data_buffer);

        // compute error from jpeg2000
        for (size_t i = 0; i < tot_size; ++i) {
            residual[i] = data[i] - decoded[i];
        }

        wavelib_forward(residual, image_dims[0], image_dims[1], WAVELET_LEVELS, &coeffs, &coeffs_size);

        if (config->residual_compression_type == QUANTILE) {
            zero_out_function(coeffs, coeffs_size, quantile);
        } else if (config->residual_compression_type == SPARSIFICATION_FACTOR) {
            double q_ratio = 1. - (1. / config->residual_cr);
            quantile = zero_out_quantile(coeffs, coeffs_size, q_ratio);
        } else if (config->residual_compression_type == MAX_ERROR) {
            double *coeffs_copy = (double *) malloc(coeffs_size * sizeof(double));
            float residual_cr = 2000;
            float cur_max_error = 0;
            do {
                memcpy(coeffs_copy, coeffs, coeffs_size * sizeof(double));
                double q_ratio = 1. - (1. / residual_cr);
                quantile = zero_out_quantile(coeffs_copy, coeffs_size, q_ratio);
                wavelib_backward(residual, image_dims[0], image_dims[1], WAVELET_LEVELS, coeffs_copy);

                cur_max_error = 0;
                for (size_t i = 0; i < tot_size; ++i) {
                    float cur_error = fabsf(data[i] - decoded[i] - residual[i]);
                    if (cur_error > cur_max_error) {
                        cur_max_error = cur_error;
                    }
                }

                residual_cr /= sqrtf(2.f);
            } while (cur_max_error > config->max_error && residual_cr >= 50);

            if (cur_max_error > config->max_error) {
                fprintf(stderr, "Could not reach error target of %f (%f instead).", config->max_error, cur_max_error);
            }
            free(coeffs);
            coeffs = coeffs_copy;
        }

        size_t zero_count = 0;
        for (size_t i = 0; i < coeffs_size; ++i) {
            if (coeffs[i] == 0) {
                zero_count++;
            }
        }

        coo_size = coeffs_size - zero_count;
        size_t coo_allocated_size = coo_size + (coeffs_size / (uint16_t)-1) + 1;
        coo_t *coo = (coo_t *) malloc(coo_allocated_size * sizeof(coo_t));

        // find min and max wavelet coefficients
        coeffs_min = coeffs_max = coeffs[0];
        for (size_t i = 1; i < coeffs_size; ++i) {
            if (coeffs[i] < coeffs_min) {
                coeffs_min = coeffs[i];
            } else if (coeffs[i] > coeffs_max) {
                coeffs_max = coeffs[i];
            }
        }

        size_t last = 0;
        size_t coo_iter = 0;
        for (size_t i = 0; i < coeffs_size; ++i) {
            double cur = coeffs[i];
            if (cur != 0 || i - last == (uint16_t)-1) {
                coo[coo_iter].c = i - last;
                assert(coo[coo_iter].c != 0 || i == 0);
                coo[coo_iter].v = ((cur - coeffs_min) / (coeffs_max - coeffs_min)) * (coo_v_t)-1;
                last = i;
                coo_iter++;
            }
        }
        assert(coo_iter <= coo_allocated_size);
        coo_size = coo_iter;

        compressed_size = ZSTD_compressBound(coo_size * sizeof(coo_t));
        compressed_coefficients = (uint8_t *) malloc(compressed_size);
        compressed_size = ZSTD_compress(compressed_coefficients, compressed_size, coo, coo_size * sizeof(coo_t), 22);

        free(coo);
        free(coeffs);
        free(residual);
        free(decoded);
    }

    size_t out_size =
            2 * sizeof(float) /* minval and maxval */ +
            sizeof(size_t) + /* coeffs size */
            2 * sizeof(double) /* coeffs_min and coeffs_max */ +
            sizeof(size_t) /* coo_size */ +
            sizeof(size_t) + compressed_size + codec_data_buffer.length;
    *out_buffer = (uint8_t *) malloc(out_size);

    uint8_t *iter = *out_buffer;
    memcpy(iter, &minval, sizeof(float));
    iter += sizeof(float);
    memcpy(iter, &maxval, sizeof(float));
    iter += sizeof(float);
    memcpy(iter, &coeffs_size, sizeof(size_t));
    iter += sizeof(size_t);
    memcpy(iter, &coeffs_min, sizeof(double));
    iter += sizeof(double);
    memcpy(iter, &coeffs_max, sizeof(double));
    iter += sizeof(double);
    memcpy(iter, &coo_size, sizeof(size_t));
    iter += sizeof(size_t);
    memcpy(iter, &compressed_size, sizeof(size_t));
    iter += sizeof(size_t);
    memcpy(iter, compressed_coefficients, compressed_size);
    iter += compressed_size;
    memcpy(iter, codec_data_buffer.buffer, codec_data_buffer.length);

    assert(iter - *out_buffer == out_size - codec_data_buffer.length);

    free(compressed_coefficients);

    codec_data_buffer_destroy(&codec_data_buffer);

    return out_size;
}

void j2k_decode_internal(float **data, size_t *height, size_t *width, float minval, float maxval, codec_data_buffer_t *codec_data_buffer) {
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
    double coeffs_min = *((double *) iter);
    iter += sizeof(double);
    double coeffs_max = *((double *) iter);
    iter += sizeof(double);
    size_t coo_size = *((size_t *) iter);
    iter += sizeof(size_t);
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

    if (compressed_coefficient_size > 0) {
        coo_t *coo = (coo_t *) malloc(coo_size * sizeof(coo_t));
        double *coeffs = (double *) calloc(coeffs_size, sizeof(double));
        float *residual = (float *) calloc(tot_size, sizeof(float));
        ZSTD_decompress(coo, coo_size * sizeof(coo_t), coefficient_data, compressed_coefficient_size);

        double *coeffs_iter = coeffs;
        for (size_t i = 0; i < coo_size; ++i) {
            coeffs_iter += coo[i].c;
            assert(coo[i].c != 0 || i == 0);
            double cur = (((double) coo[i].v / (coo_v_t)-1) * (coeffs_max - coeffs_min)) + coeffs_min;
            *coeffs_iter = cur;
        }

        wavelib_backward(residual, height, width, WAVELET_LEVELS, coeffs);

        for (size_t i = 0; i < tot_size; ++i) {
            (*out_buffer)[i] += residual[i];
        }

        free(coeffs);
        free(residual);
        free(coo);
    }

    return tot_size;
}