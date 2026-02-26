#include "hdf5_stub.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "ebcc_codec.h"
#include "log.h"

#define H5Z_FILTER_EBCC_J2K 308
#define H5Z_FILTER_EBCC_JXL 309

#if defined(EBCC_HDF5_PLUGIN_J2K) && defined(EBCC_HDF5_PLUGIN_JXL)
#error "Define exactly one of EBCC_HDF5_PLUGIN_J2K or EBCC_HDF5_PLUGIN_JXL"
#endif

#if !defined(EBCC_HDF5_PLUGIN_J2K) && !defined(EBCC_HDF5_PLUGIN_JXL)
/* Keep direct/manual compilation backward compatible by defaulting to J2K. */
#define EBCC_HDF5_PLUGIN_J2K
#endif

#if defined(EBCC_HDF5_PLUGIN_J2K)
#define H5Z_FILTER_EBCC_ID H5Z_FILTER_EBCC_J2K
#define H5Z_FILTER_EBCC_DEFAULT_COMPRESSOR BASE_COMPRESSOR_J2K
#define H5Z_FILTER_EBCC_NAME "HDF5 EBCC J2K filter"
#else
#define H5Z_FILTER_EBCC_ID H5Z_FILTER_EBCC_JXL
#define H5Z_FILTER_EBCC_DEFAULT_COMPRESSOR BASE_COMPRESSOR_JXL
#define H5Z_FILTER_EBCC_NAME "HDF5 EBCC JXL filter"
#endif

static size_t H5Z_filter_ebcc(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[],
                              size_t nbytes, size_t *buf_size, void **buf);

const H5Z_class2_t H5Z_EBCC[1] = {{
    H5Z_CLASS_T_VERS,                /* H5Z_class_t version */
    (H5Z_filter_t) H5Z_FILTER_EBCC_ID, /* Filter id number */
    1,                               /* encoder_present flag */
    1,                               /* decoder_present flag */
    H5Z_FILTER_EBCC_NAME,           /* Filter name for debugging */
    NULL,                            /* can_apply callback */
    NULL,                            /* set_local callback */
    (H5Z_func_t) H5Z_filter_ebcc,    /* Filter function */
}};

H5PL_type_t H5PLget_plugin_type(void) { return H5PL_TYPE_FILTER; }
const void *H5PLget_plugin_info(void) { return H5Z_EBCC; }

float uint_ptr_to_float(const unsigned int *ptr) {
    float value = 0.0f;
    memcpy(&value, ptr, sizeof(value));
    return value;
}

void populate_config_with_backend(codec_config_t *config,
                                  size_t cd_nelmts,
                                  const unsigned int cd_values[],
                                  size_t buf_size,
                                  base_compressor_t base_compressor) {
    config->dims[0] = buf_size / sizeof(float);
    if (config->dims[0] < cd_values[0] * cd_values[1]) {
        log_fatal("Buffer size %lu is smaller than the tile size %lu x %lu = %lu", config->dims[0], cd_values[0], cd_values[1], cd_values[0] * cd_values[1]);
        exit(1);
    }
    if (config->dims[0] % (cd_values[0] * cd_values[1]) != 0) {
        log_fatal("Buffer size %lu is not divisible by the tile size %lu x %lu = %lu", config->dims[0], cd_values[0], cd_values[1], cd_values[0] * cd_values[1]);
        exit(1);
    }
    if (cd_values[0] < 32 || cd_values[1] < 32) {
        log_fatal("Tile size %lu x %lu is too small, must be at least 32 x 32", cd_values[0], cd_values[1]);
        exit(1);
    }
    for (size_t i = 0; i < 2; ++i) {
        size_t cur = cd_values[i];
        config->dims[0] /= cur;
        config->dims[i + 1] = cur;
    }

    config->base_param = uint_ptr_to_float(&cd_values[2]);
    config->base_compressor = base_compressor;
    config->residual_compression_type = (residual_t) cd_values[3];
    switch (config->residual_compression_type) {
        case NONE:
            assert(cd_nelmts == 4);
            break;
        case MAX_ERROR:
        case RELATIVE_ERROR:
            assert(cd_nelmts == 5);
            config->error = uint_ptr_to_float(&cd_values[4]);
            break;
    }
}

void populate_config(codec_config_t *config, size_t cd_nelmts, const unsigned int cd_values[], size_t buf_size) {
    populate_config_with_backend(config, cd_nelmts, cd_values, buf_size, BASE_COMPRESSOR_J2K);
}

/**
 * cd_values:
 *  rows
 *  columns
 *  base_param as float bit-pattern
 *  residual mode: none=0, max_error=1, relative_error=2
 *  max/relative error target (as float bit-pattern) when residual mode != none
 */
static size_t H5Z_filter_ebcc(unsigned int flags,
                              size_t cd_nelmts,
                              const unsigned int cd_values[],
                              size_t nbytes,
                              size_t *buf_size,
                              void **buf) {
    if (flags & H5Z_FLAG_REVERSE) {
        float *out_buffer = NULL;
        *buf_size = ebcc_decode(*buf, nbytes, &out_buffer);

        free_buffer(*buf);
        *buf = out_buffer;

        return *buf_size;
    }

    codec_config_t config;
    populate_config_with_backend(&config, cd_nelmts, cd_values, *buf_size,
                                 (base_compressor_t) H5Z_FILTER_EBCC_DEFAULT_COMPRESSOR);

    uint8_t *out_buffer = NULL;
    *buf_size = ebcc_encode(*buf, &config, &out_buffer);

    free_buffer(*buf);
    *buf = out_buffer;

    return *buf_size;
}
