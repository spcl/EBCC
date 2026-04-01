#ifndef EBCC_CODEC_H
#define EBCC_CODEC_H

#include <stdint.h>
#include <stddef.h>

#ifndef EBCC_API
#if defined(_WIN32) && defined(EBCC_DLL_EXPORTS)
#define EBCC_API __declspec(dllexport)
#else
#define EBCC_API
#endif
#endif

#define NDIMS 3
#define EBCC_VERSION_MAJOR 0
#define EBCC_VERSION_MINOR 1

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

EBCC_API size_t ebcc_encode(float *data, codec_config_t *config, uint8_t **out_buffer);
EBCC_API size_t ebcc_decode(uint8_t *data, size_t data_size, float **out_buffer);
EBCC_API void free_buffer(void *buffer);

EBCC_API void print_config(codec_config_t *config);
EBCC_API void log_set_level_from_env(void);

#endif // EBCC_CODEC_H
