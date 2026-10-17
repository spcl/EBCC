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
#define EBCC_MIN_INTERNAL_IMAGE_DIM 32
#define EBCC_MAX_INTERNAL_IMAGE_DIM 2047
#define EBCC_VERSION_MAJOR 0
#define EBCC_VERSION_MINOR 2
/* One format for all layer combinations; version 1 remains decode-only. */
#define EBCC_HEAD_VERSION 2
#define EBCC_MAX_MCT_COMPONENTS 127 /* dense float matrix fits one MCT marker */
/* Layer selectors occupy two bits of the flags byte. Zero selects defaults. */
#define EBCC_HEADER_FLAG_CONST_FIELD 0x01
#define EBCC_HEADER_FLAG_BASE_FIXED_QUANTIZATION 0x02
#define EBCC_HEADER_FLAG_RESIDUAL_SAMPLE 0x04
#define EBCC_HEADER_FLAGS_MASK 0x07
#define EBCC_HEADER_MAGIC "EBCC"
#define EBCC_CHUNKING_HEADER_VERSION 1
#define EBCC_CHUNKING_HEADER_MAGIC "EBCK"

typedef enum {
    FIXED_CR = 0, /* Base-only encoding at the requested compression ratio. */
    MAX_ERROR,
    RELATIVE_ERROR
} residual_t;


typedef struct {
    size_t dims[NDIMS];
    float base_cr; /* Target float32 compression ratio for FIXED_CR; ignored for error bounds. */
    residual_t residual_compression_type;
    float residual_cr; /* Unused; kept for ABI compatibility. */
    float error;
    size_t chunk_dims[NDIMS];
} codec_config_t;

EBCC_API size_t ebcc_encode(float *data, codec_config_t *config, uint8_t **out_buffer);
EBCC_API size_t ebcc_decode(uint8_t *data, size_t data_size, float **out_buffer);
EBCC_API size_t ebcc_encode_chunking(float *data, codec_config_t *config, uint8_t **out_buffer);
EBCC_API size_t ebcc_encode_chunking_compat(float *data, codec_config_t *config, uint8_t **out_buffer);
EBCC_API size_t ebcc_decode_chunking(uint8_t *data, size_t data_size, float **out_buffer);
EBCC_API void free_buffer(void *buffer);

EBCC_API void print_config(codec_config_t *config);
EBCC_API void log_set_level_from_env(void);

#endif // EBCC_CODEC_H
