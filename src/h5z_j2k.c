#include <H5PLextern.h>
#include <unistd.h>
#include "j2k_codec.h"

#define H5Z_FILTER_J2K 308

static size_t H5Z_filter_j2k(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t *buf_size, void **buf);

const H5Z_class2_t H5Z_J2K[1] = {{
    H5Z_CLASS_T_VERS,               /* H5Z_class_t version */
    (H5Z_filter_t) H5Z_FILTER_J2K,  /* Filter id number             */
    1,                              /* encoder_present flag (set to true) */
    1,                              /* decoder_present flag (set to true) */
    "HDF5 j2k filter L&L",       /* Filter name for debugging    */
    NULL,                           /* The "can apply" callback     */
    NULL,                           /* The "set local" callback     */
    (H5Z_func_t) H5Z_filter_j2k,    /* The actual filter function   */
}};

H5PL_type_t H5PLget_plugin_type(void) { return H5PL_TYPE_FILTER; }
const void *H5PLget_plugin_info(void) { return H5Z_J2K; }

float uint_ptr_to_float(const unsigned int *ptr) {
    return *((float *) ptr);
}

float uint_ptr_to_double(const unsigned int *ptr) {
    return *((double *) ptr);
}

void populate_config(codec_config_t *config, size_t cd_nelmts, const unsigned int cd_values[], size_t buf_size) {
    config->dims[0] = buf_size / sizeof(float);
    for (size_t i = 0; i < 2; ++i) {
        size_t cur = cd_values[i];
        config->dims[0] /= cur;
        config->dims[i + 1] = cur;
    }

    config->base_cr = uint_ptr_to_float(&cd_values[2]);
    config->residual_compression_type = cd_values[3];
    switch (config->residual_compression_type) {
        case NONE:
            assert(cd_nelmts == 4);
            break;
        case SPARSIFICATION_FACTOR:
            assert(cd_nelmts == 5);
            config->residual_cr = uint_ptr_to_float(&cd_values[4]);
            break;
        case MAX_ERROR:
            assert(cd_nelmts == 5);
            config->max_error = uint_ptr_to_float(&cd_values[4]);
            break;
        case QUANTILE:
            assert(cd_nelmts == 6);
            config->quantile = uint_ptr_to_double(&cd_values[4]);
            break;
    }
}

/**
 * cd_values:
 *  rows
 *  columns
 *  residual mode: none=0, sparsification_ratio=1, max_error=2, quantile=3
 *  max error (as float) OR quantile (as double)
*/
// cd_values should be: frame rows, frame columns, compression_ratio, quantile (thousandth), levels, quantile (optional)
static size_t H5Z_filter_j2k(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t *buf_size, void **buf) {
    if (flags & H5Z_FLAG_REVERSE) {
        // decompress data

        float *out_buffer = NULL;
        *buf_size = decode_climate_variable(*buf, nbytes, &out_buffer);

        free(*buf);
        *buf = out_buffer;

        return *buf_size;
    } else {
        codec_config_t config;

        populate_config(&config, cd_nelmts, cd_values, *buf_size);

        uint8_t *out_buffer = NULL;
        *buf_size = encode_climate_variable(*buf, &config, &out_buffer);

        free(*buf);
        *buf = out_buffer;

        return *buf_size;
    }
}