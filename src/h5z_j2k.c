#include <H5PLextern.h>
#include <unistd.h>
#include "j2k_codec.h"

#define H5Z_FILTER_J2K 308

static size_t H5Z_filter_j2k(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t *buf_size, void **buf);

/*static htri_t can_apply_j2k(hid_t dcpl_id, hid_t type_id, hid_t space_id);*/

#ifdef J2KEMU
#define H5Z_FILTER_J2K_EMU 309
static size_t H5Z_filter_j2k_emulation(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t *buf_size, void **buf);
const H5Z_class2_t H5Z_J2K[1] = {{
    H5Z_CLASS_T_VERS,               /* H5Z_class_t version */
    (H5Z_filter_t) H5Z_FILTER_J2K_EMU,  /* Filter id number             */
    1,                              /* encoder_present flag (set to true) */
    1,                              /* decoder_present flag (set to true) */
    "HDF5 j2k filter emulation L&L",/* Filter name for debugging    */
    NULL,     /* The "can apply" callback     */
    NULL,                           /* The "set local" callback     */
    (H5Z_func_t) H5Z_filter_j2k_emulation,    /* The actual filter function   */
}};
#else
const H5Z_class2_t H5Z_J2K[1] = {{
    H5Z_CLASS_T_VERS,               /* H5Z_class_t version */
    (H5Z_filter_t) H5Z_FILTER_J2K,  /* Filter id number             */
    1,                              /* encoder_present flag (set to true) */
    1,                              /* decoder_present flag (set to true) */
    "HDF5 j2k filter L&L",          /* Filter name for debugging    */
    NULL,                           /* The "can apply" callback     */
    NULL,                           /* The "set local" callback     */
    (H5Z_func_t) H5Z_filter_j2k,    /* The actual filter function   */
}};
#endif


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

    config->base_cr = uint_ptr_to_float(&cd_values[2]);
    config->residual_compression_type = (residual_t) cd_values[3];
    switch (config->residual_compression_type) {
        case NONE:
            assert(cd_nelmts == 4);
            break;
        case SPARSIFICATION_FACTOR:
            assert(cd_nelmts == 5);
            config->residual_cr = uint_ptr_to_float(&cd_values[4]);
            break;
        case MAX_ERROR:
        case RELATIVE_ERROR:
            assert(cd_nelmts == 5);
            config->error = uint_ptr_to_float(&cd_values[4]);
            break;
        case QUANTILE:
            assert(cd_nelmts == 6);
            config->quantile = uint_ptr_to_double(&cd_values[4]);
            break;
    }
}

/* can_apply func*/
/* need hdf5_dl.c from hdf5plugins*/
/*
static htri_t can_apply_j2k(hid_t dcpl_id, hid_t type_id, hid_t space_id) {
    H5D_layout_t layout;
    H5T_class_t class;
    htri_t is_float;
    htri_t not_1d;
    htri_t is_chunked;

    layout = H5Pget_layout(dcpl_id);
    is_chunked = (layout == H5D_CHUNKED);
    class = H5Tget_class(type_id);
    is_float = (class == H5T_FLOAT);
    not_1d = (H5Sget_simple_extent_ndims(space_id) >= 2);

    return is_float && not_1d && is_chunked;
}

*/

/**
 * cd_values:
 *  rows
 *  columns
 *  residual mode: none=0, sparsification_ratio=1, max_error=2, relative_error=3, quantile=4
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

#ifdef J2KEMU
static size_t H5Z_filter_j2k_emulation(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t *buf_size, void **buf) {

    if (!(flags & H5Z_FLAG_REVERSE)) {
        codec_config_t config;

        populate_config(&config, cd_nelmts, cd_values, *buf_size);

        uint8_t *out_buffer = NULL;
        *buf_size = encode_climate_variable(*buf, &config, &out_buffer);

        free(*buf);

        float *out_buffer_f = NULL;
        *buf_size = decode_climate_variable(out_buffer, *buf_size, &out_buffer_f);
        free(out_buffer);
        *buf = out_buffer_f;
        return *buf_size;
    } else {
        // decompress data
        buf_size = 0;
        return nbytes;
    }
}
#endif