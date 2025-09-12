#ifndef EBCC_CODEC_H
#define EBCC_CODEC_H

#include <stdint.h>
#include <stddef.h>

#define NDIMS 3

typedef enum {
    NONE,
    MAX_ERROR,
    RELATIVE_ERROR
} residual_t;

typedef struct {
    size_t dims[NDIMS];
    float base_cr;
    residual_t residual_compression_type;
    float residual_cr;
    float error;
} codec_config_t;

size_t ebcc_encode(float *data, codec_config_t *config, uint8_t **out_buffer);
size_t ebcc_decode(uint8_t *data, size_t data_size, float **out_buffer);
void free_buffer(void *buffer);

void print_config(codec_config_t *config);
void log_set_level_from_env(void);

#endif // EBCC_CODEC_H