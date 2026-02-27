#ifndef EBCC_CODEC_H
#define EBCC_CODEC_H

#include <stdint.h>
#include <stddef.h>

#define NDIMS 3
#define EBCC_VERSION_MAJOR 0
#define EBCC_VERSION_MINOR 1

typedef enum {
    NONE,
    MAX_ERROR,
    RELATIVE_ERROR
} residual_t;

typedef enum {
    BASE_COMPRESSOR_J2K = 0,
    BASE_COMPRESSOR_JXL = 1
} base_compressor_t;

typedef struct {
    size_t dims[NDIMS];
    float base_param;
    base_compressor_t base_compressor;
    residual_t residual_compression_type;
    float residual_cr;
    float error;
} codec_config_t;

size_t ebcc_encode(float *data, codec_config_t *config, uint8_t **out_buffer);
size_t ebcc_decode(uint8_t *data, size_t data_size, float **out_buffer);

typedef struct {
    uint8_t *buffer;
    size_t size;        // size of the buffer (maximum allowed storage)
    size_t length;      // number of data bytes in the buffer
    size_t offset;
} codec_data_buffer_t;

void codec_data_buffer_init(codec_data_buffer_t *data);
void codec_data_buffer_rewind(codec_data_buffer_t *data);
void codec_data_buffer_clear(codec_data_buffer_t *data);
void codec_data_buffer_destroy(codec_data_buffer_t *data);

void jxl_decode_internal(float **data,
                         size_t *height,
                         size_t *width,
                         float minval,
                         float maxval,
                         codec_data_buffer_t *codec_data_buffer);
size_t jxl_encode_internal_with_effort(void *data,
                                       size_t *image_dims,
                                       float distance,
                                       codec_data_buffer_t *codec_data_buffer,
                                       int effort);
void free_buffer(void *buffer);

void print_config(codec_config_t *config);
void log_set_level_from_env(void);

#endif // EBCC_CODEC_H
