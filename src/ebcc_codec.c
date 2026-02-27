#include "ebcc_codec.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "log.h"
#ifdef ENABLE_PERF
#include <sys/prctl.h>
#endif

#include "openjpeg.h"
#include "zstd.h"
#include <assert.h>
#ifdef ENABLE_JXL
#include <jxl/encode.h>
#include <jxl/decode.h>
#endif
#ifdef __linux__
#include <malloc.h>
#else
#include <stdlib.h>
#endif
#include "spiht_re.h"


#define MIN(x, y) ((x)<(y)?(x):(y))
#define CEIL(x,y) (((x) + (y) - 1) / (y))
#define TRUE 1
#define FALSE 0

#define WAVELET_LEVELS 3
#define EBCC_HEADER_VERSION 2
#define EBCC_HEADER_FLAG_CONST_FIELD 0x01
#define EBCC_HEADER_MAGIC "EBCC"
#define J2K_PARAM_MIN 1.0f
#define J2K_PARAM_MAX 1000.0f
#define JXL_PARAM_MIN 0.0f
#define JXL_PARAM_MAX 25.0f

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
#define EBCC_STATIC_ASSERT(cond, msg) _Static_assert(cond, msg)
#else
#define EBCC_STATIC_ASSERT_CAT_(a, b) a##b
#define EBCC_STATIC_ASSERT_CAT(a, b) EBCC_STATIC_ASSERT_CAT_(a, b)
#define EBCC_STATIC_ASSERT(cond, msg) typedef char EBCC_STATIC_ASSERT_CAT(ebcc_static_assert_, __LINE__)[(cond) ? 1 : -1]
#endif

typedef struct {
    size_t (*encode)(void *data, size_t *image_dims, size_t *tile_dims, float param, codec_data_buffer_t *out);
    void (*decode)(float **data, size_t *height, size_t *width, float minval, float maxval, codec_data_buffer_t *stream);
    float (*error_bound_search)(uint16_t *scaled_data, size_t *image_dims, size_t *tile_dims, float param,
                                codec_data_buffer_t *codec_data_buffer, float **decoded, float minval, float maxval,
                                float *data, size_t tot_size, float error_target, double base_quantile_target);
    const char *name;
    uint8_t id;
    const char *param_name;
    float param_min;
    float param_max;
} base_codec_backend_t;

static size_t j2k_encode_internal(void *data, size_t *image_dims, size_t *tile_dims, float base_param,
                                  codec_data_buffer_t *codec_data_buffer);
static void j2k_decode_internal(float **data, size_t *height, size_t *width, float minval, float maxval,
                                codec_data_buffer_t *codec_data_buffer);
static size_t jxl_encode_internal(void *data, size_t *image_dims, size_t *tile_dims, float base_param,
                                  codec_data_buffer_t *codec_data_buffer);
void jxl_decode_internal(float **data, size_t *height, size_t *width, float minval, float maxval,
                                codec_data_buffer_t *codec_data_buffer);

void codec_data_buffer_init(codec_data_buffer_t* data) {
    const size_t initial_buffer_size = 1024;
    data->buffer = (uint8_t *) malloc(initial_buffer_size);
    data->size = initial_buffer_size;
    data->length = 0;
    data->offset = 0;
}

void codec_data_buffer_rewind(codec_data_buffer_t* data) {
    data->offset = 0;
}

void codec_data_buffer_clear(codec_data_buffer_t* data) {
    data->length = 0;
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

static size_t write_to_buffer_stream(void *input_buffer, size_t len, codec_data_buffer_t *stream_data) {
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

static size_t read_from_buffer_stream(void *output_buffer, size_t len, codec_data_buffer_t *stream_data) {
    if (stream_data->offset >= stream_data->length) {
        return (size_t) -1;
    }

    size_t n_bytes_to_read = MIN(len, stream_data->length - stream_data->offset);
    memcpy(output_buffer, stream_data->buffer + stream_data->offset, n_bytes_to_read);
    stream_data->offset += n_bytes_to_read;

    return n_bytes_to_read;
}

static size_t j2k_encode_internal(void *data, size_t *image_dims, size_t *tile_dims, float base_param,
        codec_data_buffer_t *codec_data_buffer) {
    size_t n_tiles = image_dims[0] / tile_dims[0];
    size_t tile_size = tile_dims[0] * tile_dims[1];

    opj_cparameters_t parameters;
    opj_set_default_encoder_parameters(&parameters);

    // Set image parameters
    parameters.tcp_numlayers = 1;
    parameters.cp_disto_alloc = 1;
    parameters.tcp_rates[0] = base_param / 2;
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
    return codec_data_buffer->length;
}

typedef uint8_t coo_v_t;

typedef struct {
    uint16_t c;
    coo_v_t v;
} coo_t;


typedef struct {
    uint8_t magic[4];
    uint8_t version;
    uint8_t flags;
    uint8_t base_compressor;
    uint8_t reserved;
    uint32_t minval_bits;
    uint32_t maxval_bits;
    uint64_t coeffs_size;
    uint32_t residual_minval_bits;
    uint32_t residual_maxval_bits;
    uint64_t compressed_size;
    uint64_t tail_size;
} ebcc_header_t;

EBCC_STATIC_ASSERT(sizeof(float) == 4, "EBCC serialization requires 32-bit float");
EBCC_STATIC_ASSERT(sizeof(ebcc_header_t) == 48, "EBCC header must be fixed-size");

static uint32_t float_to_bits(float value) {
    uint32_t bits = 0;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

static float bits_to_float(uint32_t bits) {
    float value = 0.0f;
    memcpy(&value, &bits, sizeof(value));
    return value;
}

static void assert_endianness(void) {
    uint32_t value = 0x01020304u;
    uint8_t bytes[sizeof(uint32_t)];
    memcpy(bytes, &value, sizeof(bytes));
    assert(bytes[0] == 0x04);
}

const char* residual_t_names[] = {
    "NONE",
    "MAX_ERROR",
    "RELATIVE_ERROR"
};

static const char* base_compressor_names[] = {
    "j2k",
    "jxl"
};

static const char *base_compressor_name(base_compressor_t compressor) {
    if (compressor < BASE_COMPRESSOR_J2K || compressor > BASE_COMPRESSOR_JXL) {
        return "unknown";
    }
    return base_compressor_names[(int) compressor];
}

void print_config(codec_config_t *config) {
    log_info("dimensions:\t(%lu, %lu, %lu)", config->dims[0], config->dims[1], config->dims[2]);
    log_info("base_param:\t%f", config->base_param);
    log_info("base compressor:\t%s", base_compressor_name(config->base_compressor));
    log_info("residual type:\t%s", residual_t_names[config->residual_compression_type]);
    switch (config->residual_compression_type) {
        case NONE:
            break;
        case MAX_ERROR:
            log_info("max error:\t%f", config->error);
            break;
        case RELATIVE_ERROR:
            log_info("relative error:\t%f", config->error);
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

double get_mean_error(const float *data, const float *decoded, const float *residual, const size_t tot_size) {
    double sum_error = 0;
    for (size_t i = 0; i < tot_size; ++i) {
        float residual_value = residual ? residual[i] : 0;
        sum_error += data[i] - (decoded[i] + residual_value);
    }
    return (sum_error / tot_size);
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

static int get_env_int(const char *name, int default_value, int min_value, int max_value) {
    const char *raw = getenv(name);
    if (!raw) {
        return default_value;
    }
    char *endptr = NULL;
    long parsed = strtol(raw, &endptr, 10);
    if (*endptr != '\0' || parsed < min_value || parsed > max_value) {
        log_warn("Ignoring invalid %s=%s; expected integer in [%d, %d]", name, raw, min_value, max_value);
        return default_value;
    }
    return (int) parsed;
}

size_t jxl_encode_internal_with_effort(void *data,
                                              size_t *image_dims,
                                              float distance,
                                              codec_data_buffer_t *codec_data_buffer,
                                              int effort) {
#ifdef ENABLE_JXL
    if (distance < JXL_PARAM_MIN) {
        distance = JXL_PARAM_MIN;
    } else if (distance > JXL_PARAM_MAX) {
        distance = JXL_PARAM_MAX;
    }

    JxlEncoder *encoder = JxlEncoderCreate(NULL);
    if (!encoder) {
        log_fatal("Failed to create JXL encoder");
        return 0;
    }

    JxlEncoderUseContainer(encoder, JXL_FALSE);

    JxlBasicInfo basic_info;
    JxlEncoderInitBasicInfo(&basic_info);
    basic_info.xsize = (uint32_t) image_dims[1];
    basic_info.ysize = (uint32_t) image_dims[0];
    basic_info.bits_per_sample = 16;
    basic_info.exponent_bits_per_sample = 0;
    basic_info.num_color_channels = 1;
    basic_info.num_extra_channels = 0;
    basic_info.uses_original_profile = distance == 0.0f ? JXL_TRUE : JXL_FALSE; // set to true affects compression ratio !
    if (JxlEncoderSetBasicInfo(encoder, &basic_info) != JXL_ENC_SUCCESS) {
        log_fatal("Failed to configure JXL basic info");
        JxlEncoderDestroy(encoder);
        return 0;
    }

    JxlColorEncoding color_encoding;
    JxlColorEncodingSetToLinearSRGB(&color_encoding, JXL_TRUE);
    if (JxlEncoderSetColorEncoding(encoder, &color_encoding) != JXL_ENC_SUCCESS) {
        log_fatal("Failed to configure JXL color encoding");
        JxlEncoderDestroy(encoder);
        return 0;
    }

    JxlEncoderFrameSettings *frame_settings = JxlEncoderFrameSettingsCreate(encoder, NULL);
    if (!frame_settings) {
        log_fatal("Failed to create JXL frame settings");
        JxlEncoderDestroy(encoder);
        return 0;
    }
    JxlEncoderFrameSettingsSetOption(frame_settings, JXL_ENC_FRAME_SETTING_EFFORT, effort);
    JxlEncoderSetFrameDistance(frame_settings, distance);
    if (distance == 0.0f &&
        JXL_ENC_SUCCESS != JxlEncoderSetFrameLossless(frame_settings, JXL_TRUE)) {
      log_fatal("JxlEncoderSetFrameLossless() failed.");
      JxlEncoderDestroy(encoder);
      return 0;
    }

    JxlPixelFormat pixel_format = {1, JXL_TYPE_UINT16, JXL_NATIVE_ENDIAN, 0};
    size_t image_size_bytes = image_dims[0] * image_dims[1] * sizeof(uint16_t);
    if (JxlEncoderAddImageFrame(frame_settings, &pixel_format, data, image_size_bytes) != JXL_ENC_SUCCESS) {
        log_fatal("Failed to add JXL frame");
        JxlEncoderDestroy(encoder);
        return 0;
    }
    JxlEncoderCloseInput(encoder);

    codec_data_buffer_clear(codec_data_buffer);

    uint8_t jxl_chunk[65536];
    uint8_t *next_out = jxl_chunk;
    size_t avail_out = sizeof(jxl_chunk);
    for (;;) {
        JxlEncoderStatus status = JxlEncoderProcessOutput(encoder, &next_out, &avail_out);
        size_t produced = sizeof(jxl_chunk) - avail_out;
        if (produced > 0) {
            write_to_buffer_stream(jxl_chunk, produced, codec_data_buffer);
        }
        if (status == JXL_ENC_SUCCESS) {
            break;
        }
        if (status != JXL_ENC_NEED_MORE_OUTPUT) {
            log_fatal("JXL encoder failed while producing output");
            JxlEncoderDestroy(encoder);
            return 0;
        }
        next_out = jxl_chunk;
        avail_out = sizeof(jxl_chunk);
    }
    JxlEncoderDestroy(encoder);
    return codec_data_buffer->length;
#else
    (void) data;
    (void) image_dims;
    (void) distance;
    (void) codec_data_buffer;
    (void) effort;
    log_fatal("JXL support not compiled. Rebuild with -DENABLE_JXL=ON.");
    return 0;
#endif
}

static size_t jxl_encode_internal(void *data,
                                  size_t *image_dims,
                                  size_t *tile_dims,
                                  float base_param,
                                  codec_data_buffer_t *codec_data_buffer) {
    (void) tile_dims;
    int effort = get_env_int("EBCC_JXL_EFFORT", 7, 1, 9);
    return jxl_encode_internal_with_effort(data, image_dims, base_param, codec_data_buffer, effort);
}

void jxl_decode_internal(float **data,
                                size_t *height,
                                size_t *width,
                                float minval,
                                float maxval,
                                codec_data_buffer_t *codec_data_buffer) {
#ifdef ENABLE_JXL
    codec_data_buffer_rewind(codec_data_buffer);

    JxlDecoder *decoder = JxlDecoderCreate(NULL);
    if (!decoder) {
        log_fatal("Failed to create JXL decoder");
        return;
    }
    if (JxlDecoderSubscribeEvents(decoder, JXL_DEC_BASIC_INFO | JXL_DEC_FULL_IMAGE) != JXL_DEC_SUCCESS) {
        log_fatal("Failed to subscribe JXL decoder events");
        JxlDecoderDestroy(decoder);
        return;
    }

    JxlDecoderSetInput(decoder, codec_data_buffer->buffer, codec_data_buffer->length);
    JxlDecoderCloseInput(decoder);

    JxlPixelFormat pixel_format = {1, JXL_TYPE_UINT16, JXL_NATIVE_ENDIAN, 0};
    uint16_t *decoded_u16 = NULL;
    size_t decoded_u16_size = 0;
    size_t xsize = 0, ysize = 0;

    for (;;) {
        JxlDecoderStatus status = JxlDecoderProcessInput(decoder);
        if (status == JXL_DEC_ERROR) {
            log_fatal("JXL decoder failed");
            free(decoded_u16);
            JxlDecoderDestroy(decoder);
            return;
        }
        if (status == JXL_DEC_NEED_MORE_INPUT) {
            log_fatal("Truncated JXL payload");
            free(decoded_u16);
            JxlDecoderDestroy(decoder);
            return;
        }
        if (status == JXL_DEC_BASIC_INFO) {
            JxlBasicInfo info;
            if (JxlDecoderGetBasicInfo(decoder, &info) != JXL_DEC_SUCCESS) {
                log_fatal("Failed to read JXL basic info");
                free(decoded_u16);
                JxlDecoderDestroy(decoder);
                return;
            }
            xsize = info.xsize;
            ysize = info.ysize;
            if (width) {
                *width = xsize;
            }
            if (height) {
                *height = ysize;
            }
            if (JxlDecoderImageOutBufferSize(decoder, &pixel_format, &decoded_u16_size) != JXL_DEC_SUCCESS) {
                log_fatal("Failed to query JXL output buffer size");
                free(decoded_u16);
                JxlDecoderDestroy(decoder);
                return;
            }
            decoded_u16 = (uint16_t *) malloc(decoded_u16_size);
            if (!decoded_u16) {
                log_fatal("Failed to allocate JXL output buffer");
                JxlDecoderDestroy(decoder);
                return;
            }
            if (JxlDecoderSetImageOutBuffer(decoder, &pixel_format, decoded_u16, decoded_u16_size) != JXL_DEC_SUCCESS) {
                log_fatal("Failed to set JXL output buffer");
                free(decoded_u16);
                JxlDecoderDestroy(decoder);
                return;
            }
            continue;
        }
        if (status == JXL_DEC_FULL_IMAGE) {
            continue;
        }
        if (status == JXL_DEC_SUCCESS) {
            break;
        }
    }

    size_t num_pixels = xsize * ysize;
    if (!*data) {
        *data = (float *) malloc(num_pixels * sizeof(float));
    }
    for (size_t i = 0; i < num_pixels; ++i) {
        (*data)[i] = ((float) decoded_u16[i] / (uint16_t)-1) * (maxval - minval) + minval;
    }

    free(decoded_u16);
    JxlDecoderDestroy(decoder);
#else
    (void) data;
    (void) height;
    (void) width;
    (void) minval;
    (void) maxval;
    (void) codec_data_buffer;
    log_fatal("JXL support not compiled. Rebuild with -DENABLE_JXL=ON.");
#endif
}

double emulate_j2k_compression(uint16_t *scaled_data, size_t *image_dims, size_t *tile_dims, float current_cr,
                             codec_data_buffer_t *codec_data_buffer, float **decoded, float minval, float maxval,
                             float *data, size_t tot_size, float error_target) {
    codec_data_buffer_clear(codec_data_buffer);
    j2k_encode_internal(scaled_data, image_dims, tile_dims, current_cr, codec_data_buffer);
    codec_data_buffer_rewind(codec_data_buffer);
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

double emulate_jxl_compression(uint16_t *scaled_data, size_t *image_dims, size_t *tile_dims, float current_distance,
                             codec_data_buffer_t *codec_data_buffer, float **decoded, float minval, float maxval,
                             float *data, size_t tot_size, float error_target) {
    (void) tile_dims;
    int effort = get_env_int("EBCC_JXL_SEARCH_EFFORT", 4, 1, 9);
    codec_data_buffer_clear(codec_data_buffer);
    jxl_encode_internal_with_effort(scaled_data, image_dims, current_distance, codec_data_buffer, effort);
    codec_data_buffer_rewind(codec_data_buffer);
    jxl_decode_internal(decoded, NULL, NULL, minval, maxval, codec_data_buffer);
    return get_error_target_quantile(data, *decoded, NULL, tot_size, error_target);
}

float error_bound_jxl_compression(uint16_t *scaled_data, size_t *image_dims, size_t *tile_dims, float current_distance,
                             codec_data_buffer_t *codec_data_buffer, float **decoded, float minval, float maxval,
                             float *data, size_t tot_size, float error_target, double base_quantile_target) {
    float dist_lo = JXL_PARAM_MIN;
    float dist_hi = JXL_PARAM_MAX;
    float best_feasible_dist = dist_lo;
    double eps = 1e-8;

    double quantile_at_lo = emulate_jxl_compression(scaled_data, image_dims, tile_dims, dist_lo,
                                                    codec_data_buffer, decoded, minval, maxval,
                                                    data, tot_size, error_target);
    log_trace("current_distance: %f, 1-error_target_quantile: %.1e, jxl_length: %lu",
              current_distance, 1 - quantile_at_lo, codec_data_buffer->length);
    log_trace("dist_lo: %f, 1-error_target_quantile: %.1e, jxl_length: %lu",
              dist_lo, 1 - quantile_at_lo, codec_data_buffer->length);
    if (quantile_at_lo < base_quantile_target) {
        log_warn("JXL distance %.3f does not satisfy the error target quantile; using minimum distance.", dist_lo);
        return dist_lo;
    }

    double quantile_at_hi = emulate_jxl_compression(scaled_data, image_dims, tile_dims, dist_hi,
                                                    codec_data_buffer, decoded, minval, maxval,
                                                    data, tot_size, error_target);
    log_trace("dist_hi: %f, 1-error_target_quantile: %.1e, jxl_length: %lu",
              dist_hi, 1 - quantile_at_hi, codec_data_buffer->length);
    if (quantile_at_hi >= base_quantile_target) {
        return dist_hi;
    }

    while (dist_hi - dist_lo > 1e-3f) {
        float mid = 0.5f * (dist_hi + dist_lo);
        double quantile_mid = emulate_jxl_compression(scaled_data, image_dims, tile_dims, mid,
                                                      codec_data_buffer, decoded, minval, maxval,
                                                      data, tot_size, error_target);
        log_trace("current_distance: %f, dist_lo: %f, dist_hi: %f, 1-error_target_quantile: %.1e, jxl_length: %lu",
                  mid, dist_lo, dist_hi, 1 - quantile_mid, codec_data_buffer->length);
        if (quantile_mid + eps >= base_quantile_target) {
            best_feasible_dist = mid;
            dist_lo = mid;
        } else {
            dist_hi = mid;
        }
    }
    {
        double best_quantile = emulate_jxl_compression(scaled_data, image_dims, tile_dims, best_feasible_dist,
                                                       codec_data_buffer, decoded, minval, maxval,
                                                       data, tot_size, error_target);
        log_trace("best_distance: %f, 1-error_target_quantile: %.1e, jxl_length: %lu",
                  best_feasible_dist, 1 - best_quantile, codec_data_buffer->length);
    }
    return best_feasible_dist;
}

static const base_codec_backend_t j2k_backend = {
    .encode = j2k_encode_internal,
    .decode = j2k_decode_internal,
    .error_bound_search = error_bound_j2k_compression,
    .name = "j2k",
    .id = BASE_COMPRESSOR_J2K,
    .param_name = "cr",
    .param_min = J2K_PARAM_MIN,
    .param_max = J2K_PARAM_MAX,
};

static const base_codec_backend_t jxl_backend = {
    .encode = jxl_encode_internal,
    .decode = jxl_decode_internal,
    .error_bound_search = error_bound_jxl_compression,
    .name = "jxl",
    .id = BASE_COMPRESSOR_JXL,
    .param_name = "distance",
    .param_min = JXL_PARAM_MIN,
    .param_max = JXL_PARAM_MAX,
};

static const base_codec_backend_t *get_backend_by_id(uint8_t id) {
    if (id == BASE_COMPRESSOR_J2K) {
        return &j2k_backend;
    }
    if (id == BASE_COMPRESSOR_JXL) {
        return &jxl_backend;
    }
    return NULL;
}

static const base_codec_backend_t *get_backend_for_encode(codec_config_t *config) {
#ifndef ENABLE_JXL
    if (config->base_compressor == BASE_COMPRESSOR_JXL) {
        log_fatal("JXL support not compiled. Rebuild with -DENABLE_JXL=ON.");
        return NULL;
    }
#endif
    const base_codec_backend_t *backend = get_backend_by_id((uint8_t) config->base_compressor);
    if (!backend) {
        log_fatal("Unknown base compressor: %d", (int) config->base_compressor);
        return NULL;
    }
    if (config->base_param < backend->param_min || config->base_param > backend->param_max) {
        log_warn("Base %s parameter %f is outside [%f, %f]", backend->param_name, config->base_param,
                 backend->param_min, backend->param_max);
    }
    return backend;
}

void check_nan_inf(const float *data, const size_t tot_size) {
    for (size_t i = 0; i < tot_size; ++i) {
        if (isnan(data[i]) || isinf(data[i])) {
            log_fatal("NaN or Inf found in data at index %lu", i);
            exit(1);
        }
    }
}

size_t ebcc_encode(float *data, codec_config_t *config, uint8_t **out_buffer) {
#ifdef ENABLE_PERF
    prctl(PR_TASK_PERF_EVENTS_ENABLE);
#endif
    assert_endianness();


    int pure_base_codec_required = FALSE, pure_base_codec_done = FALSE, pure_base_codec_disabled = FALSE, pure_base_codec_consistency_disabled = FALSE, mean_error_adjustment_disabled = FALSE;
    size_t compressed_size = 0, base_codec_buffer_length = 0;
    uint8_t *compressed_coefficients = NULL;
    uint8_t *coeffs_buf = NULL;
    uint8_t *base_codec_buffer = NULL; 
    float residual_maxval = 0., residual_minval = 0., error_target = -1, current_base_param = -1;
    size_t coeffs_size = 0, coeffs_size_orig = 0, coeffs_trunc_bits = 0; /*coeffs_size: #bytes*/
    double trunc_hi, trunc_lo, best_feasible_trunc;
    double eps = 1e-8, base_error_quantile=1e-6;
    double cur_mean_error = 0.0;

    // Load base_error_quantile from env var EBCC_INIT_BASE_ERROR_QUANTILE, default to 1e-6 if not set

    log_set_level_from_env();
    
    const char *env_base_error_quantile = getenv("EBCC_INIT_BASE_ERROR_QUANTILE");
    const char *env_disable_pure_base_codec_fallback = getenv("EBCC_DISABLE_PURE_BASE_COMPRESSION_FALLBACK");
    const char *env_disable_pure_base_codec_consistency = getenv("EBCC_DISABLE_PURE_BASE_COMPRESSION_FALLBACK_CONSISTENCY");
    const char *env_disable_mean_adjustment = getenv("EBCC_DISABLE_MEAN_ADJUSTMENT");
    if (env_base_error_quantile) {
        base_error_quantile = strtod(env_base_error_quantile, NULL);
    }
    if (env_disable_pure_base_codec_fallback) {
        pure_base_codec_disabled = TRUE;
    }
    if (env_disable_pure_base_codec_consistency) {
        pure_base_codec_consistency_disabled = TRUE;
    }
    if (env_disable_mean_adjustment) {
        mean_error_adjustment_disabled = TRUE;
    }
    double base_quantile_target = 1 - base_error_quantile;

    print_config(config);
    log_info("1 - base_quantile_target: %.1e", 1-base_quantile_target);
    log_info("Disable pure base compression fallback: %s", pure_base_codec_disabled ? "TRUE" : "FALSE");

    const base_codec_backend_t *backend = get_backend_for_encode(config);
    if (!backend) {
        return 0;
    }
    log_info("[ebcc] using backend: %s", backend->name);



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

        // encode using selected backend
        codec_data_buffer_clear(&codec_data_buffer);
        if (backend->encode(scaled_data, image_dims, tile_dims, config->base_param, &codec_data_buffer) == 0) {
            free(scaled_data);
            codec_data_buffer_destroy(&codec_data_buffer);
            return 0;
        }
        if (config->residual_compression_type == NONE) {
            base_codec_buffer = (uint8_t *) malloc(codec_data_buffer.length);
            memcpy(base_codec_buffer, codec_data_buffer.buffer, codec_data_buffer.length);
            base_codec_buffer_length = codec_data_buffer.length;
        }
        codec_data_buffer_rewind(&codec_data_buffer);
    }

    if ((config->residual_compression_type != NONE) && !const_field) {
        // decode back the image
        float *residual = (float *) malloc(tot_size * sizeof(float));
        float *residual_norm = (float *) malloc(tot_size * sizeof(float));
        float *decoded = (float *) malloc(tot_size * sizeof(float));
        backend->decode(&decoded, NULL, NULL, minval, maxval, &codec_data_buffer);

        cur_mean_error = get_mean_error(data, decoded, NULL, tot_size);

        // compute error from jpeg2000
        for (size_t i = 0; i < tot_size; ++i) {
            residual[i] = data[i] - decoded[i];
        }

        findMinMaxf(residual, tot_size, &residual_minval, &residual_maxval);
        for (size_t i = 0; i < tot_size; ++i) {
            residual_norm[i] = (residual[i] - residual_minval) / (residual_maxval - residual_minval);
        }

        if (config->residual_compression_type == MAX_ERROR ||
                config->residual_compression_type == RELATIVE_ERROR) {
            error_target = config->error, current_base_param = config->base_param;
            if (config->residual_compression_type == RELATIVE_ERROR) {
                error_target *= get_data_range(data, tot_size);
            }

            current_base_param = backend->error_bound_search(scaled_data, image_dims, tile_dims, current_base_param,
                                                             &codec_data_buffer, &decoded, minval, maxval,
                                                             data, tot_size, error_target, base_quantile_target);
            
            for (size_t i = 0; i < tot_size; ++i) {
                residual[i] = data[i] - decoded[i];
            }
            findMinMaxf(residual, tot_size, &residual_minval, &residual_maxval);

            float cur_max_error = fmaxf(fabsf(residual_minval), fabsf(residual_maxval));
            float best_feasible_error = -1;
            int skip_residual = cur_max_error <= error_target;
            pure_base_codec_done = base_quantile_target == 1.0;

            if (pure_base_codec_done) log_info("Pure base compression is feasible, compression error: %f, %s: %f", cur_max_error, backend->param_name, current_base_param);
            

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
                    log_info("Could not reach error target of %f (%f instead), base compression max error: %f. Retry with pure base compression.", error_target, cur_max_error, fmaxf(fabsf(residual_minval), fabsf(residual_maxval)));
                    skip_residual = TRUE;
                    pure_base_codec_required = TRUE;
                    /*DONE: if this happen, go for full jpeg2000*/
                } else {
                    best_feasible_error = cur_max_error;
                    cur_mean_error = get_mean_error(data, decoded, residual, tot_size);
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
                            cur_mean_error = get_mean_error(data, decoded, residual, tot_size);
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

        /* Try again with pure base compression, to see if adding residual compression has higher compression ratio */
        base_codec_buffer_length = codec_data_buffer.length;
        size_t base_codec_buffer_size_limit = 2 * (compressed_size + base_codec_buffer_length);
        base_codec_buffer = (uint8_t *) malloc(base_codec_buffer_size_limit);
        memcpy(base_codec_buffer, codec_data_buffer.buffer, base_codec_buffer_length);
        if ((! pure_base_codec_done) && (! pure_base_codec_disabled) && (config->residual_compression_type == MAX_ERROR ||
            config->residual_compression_type == RELATIVE_ERROR)) {
            assert(error_target > 0);
            /* ===========Maintain consistency with quantile = 0 (Not necessary) =========== */
            if (!pure_base_codec_consistency_disabled) {
                codec_data_buffer_clear(&codec_data_buffer);
                if (backend->encode(scaled_data, image_dims, tile_dims, config->base_param, &codec_data_buffer) == 0) {
                    free(coeffs_buf);
                    free(residual);
                    free(residual_norm);
                    free(decoded);
                    free(scaled_data);
                    codec_data_buffer_destroy(&codec_data_buffer);
                    return 0;
                }
                codec_data_buffer_rewind(&codec_data_buffer);
                backend->decode(&decoded, NULL, NULL, minval, maxval, &codec_data_buffer);
                current_base_param = config->base_param;
            }
            /* ===========Maintain consistency with quantile = 0 (Not necessary) =========== */
            backend->error_bound_search(scaled_data, image_dims, tile_dims, current_base_param,
                                        &codec_data_buffer, &decoded, minval, maxval,
                                        data, tot_size, error_target, 1.0);

            if ((codec_data_buffer.length < compressed_size + base_codec_buffer_length) || pure_base_codec_required) {
                /* Pure JP2 is better than JP2 + SPWV */
                if (codec_data_buffer.length < compressed_size + base_codec_buffer_length)
                    log_info("Pure base compression (%lu) is better than base (%lu) + residual (%lu) compression (sum: %lu)", codec_data_buffer.length, base_codec_buffer_length, compressed_size, compressed_size + base_codec_buffer_length);
                
                cur_mean_error = get_mean_error(data, decoded, NULL, tot_size);

                compressed_size = 0;
                coeffs_size = 0;
                if (codec_data_buffer.length > base_codec_buffer_size_limit) { /* This can happen when pure_base_codec_required enabled */
                    free(base_codec_buffer);
                    base_codec_buffer = (uint8_t *) malloc(codec_data_buffer.length);
                }
                base_codec_buffer_length = codec_data_buffer.length;
                memcpy(base_codec_buffer, codec_data_buffer.buffer, base_codec_buffer_length);
            }
        }
        free(coeffs_buf);
        free(residual);
        free(residual_norm);
        free(decoded);
    }

    if (!const_field) free(scaled_data);

    log_info("Mean of compression error: %e", cur_mean_error);
    if (!mean_error_adjustment_disabled && fabs(cur_mean_error) > 1e-18) {
        minval += cur_mean_error;
        maxval += cur_mean_error;
        log_info("Adjusting minval and maxval to %f, %f", minval, maxval);
    }

    size_t codec_size = (const_field) ? sizeof(uint64_t) : base_codec_buffer_length; /*Only output array length if having a constant field*/

    size_t out_size =
            sizeof(ebcc_header_t) +
            compressed_size + codec_size;
    *out_buffer = (uint8_t *) malloc(out_size);

    log_info("coeffs_size: %lu, compressed_size: %lu, base_length: %lu, compression ratio: %f", coeffs_size, compressed_size, codec_size, (double) (tot_size * sizeof(float)) / out_size);

    uint8_t *iter = *out_buffer;
    ebcc_header_t header = {0};
    memcpy(header.magic, EBCC_HEADER_MAGIC, 4);
    header.version = EBCC_HEADER_VERSION;
    if (const_field) {
        header.flags |= EBCC_HEADER_FLAG_CONST_FIELD;
    }
    header.base_compressor = backend->id;
    header.minval_bits = float_to_bits(minval);
    header.maxval_bits = float_to_bits(maxval);
    header.coeffs_size = (uint64_t) coeffs_size;
    header.residual_minval_bits = float_to_bits(residual_minval);
    header.residual_maxval_bits = float_to_bits(residual_maxval);
    header.compressed_size = (uint64_t) compressed_size;
    header.tail_size = (uint64_t) codec_size;
    memcpy(iter, &header, sizeof(header));
    iter += sizeof(header);
    if (compressed_size > 0) {
        memcpy(iter, compressed_coefficients, compressed_size);
        iter += compressed_size;
    }
    if (const_field) {
        uint64_t tot_size_u64 = (uint64_t) tot_size;
        memcpy(iter, &tot_size_u64, sizeof(uint64_t));
        iter += sizeof(uint64_t);
    } else {
        memcpy(iter, base_codec_buffer, base_codec_buffer_length);
        iter += base_codec_buffer_length;
    }
    assert((size_t) (iter - *out_buffer) == out_size);

    if (compressed_coefficients) free(compressed_coefficients);
    if (base_codec_buffer) free(base_codec_buffer);

    codec_data_buffer_destroy(&codec_data_buffer);

#ifdef ENABLE_PERF
    prctl(PR_TASK_PERF_EVENTS_DISABLE);
#endif
    return out_size;
}

static void j2k_decode_internal(float **data, size_t *height, size_t *width, float minval, float maxval,
        codec_data_buffer_t *codec_data_buffer) {
    codec_data_buffer_rewind(codec_data_buffer);

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

static int legacy_read_bytes(const uint8_t **iter, const uint8_t *end, void *dst, size_t nbytes) {
    if ((size_t) (end - *iter) < nbytes) {
        return FALSE;
    }
    memcpy(dst, *iter, nbytes);
    *iter += nbytes;
    return TRUE;
}

size_t ebcc_decode_legacy(uint8_t *data, size_t data_size, float **out_buffer) {
    codec_data_buffer_t codec_data_buffer;
    size_t tot_size = 0;
    const uint8_t *iter = data;
    const uint8_t *end = data + data_size;

    float minval = 0.0f, maxval = 0.0f, residual_minval = 0.0f, residual_maxval = 0.0f;
    size_t coeffs_size = 0, compressed_coefficient_size = 0;
    if (!legacy_read_bytes(&iter, end, &minval, sizeof(float)) ||
            !legacy_read_bytes(&iter, end, &maxval, sizeof(float)) ||
            !legacy_read_bytes(&iter, end, &coeffs_size, sizeof(size_t)) ||
            !legacy_read_bytes(&iter, end, &residual_minval, sizeof(float)) ||
            !legacy_read_bytes(&iter, end, &residual_maxval, sizeof(float)) ||
            !legacy_read_bytes(&iter, end, &compressed_coefficient_size, sizeof(size_t))) {
        log_fatal("Invalid legacy encoded data: truncated header");
        return 0;
    }

    if ((size_t) (end - iter) < compressed_coefficient_size) {
        log_fatal("Invalid legacy encoded data: truncated residual payload");
        return 0;
    }
    uint8_t *coefficient_data = (uint8_t *) iter;
    iter += compressed_coefficient_size;

    size_t height = 0, width = 0;
    int const_field = minval == maxval;
    if (const_field) {
        if (!legacy_read_bytes(&iter, end, &tot_size, sizeof(size_t))) {
            log_fatal("Invalid legacy encoded data: missing const-field length");
            return 0;
        }
        *out_buffer = (float *) malloc(tot_size * sizeof(float));
        for (size_t i = 0; i < tot_size; ++i) {
            (*out_buffer)[i] = minval;
        }
    } else {
        codec_data_buffer.buffer = (uint8_t *) iter;
        codec_data_buffer.size = (size_t) (end - iter);
        codec_data_buffer.length = (size_t) (end - iter);
        codec_data_buffer.offset = 0;
        j2k_decode_internal(out_buffer, &height, &width, minval, maxval, &codec_data_buffer);
        tot_size = height * width;
    }

    log_info("Compression Ratio: %f", (double) (tot_size * sizeof(float)) / data_size);

    if (compressed_coefficient_size > 0 && coeffs_size > 0) {
        if (const_field) {
            log_fatal("Invalid legacy encoded data: residual data cannot be applied to const field");
            return 0;
        }
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

size_t ebcc_decode(uint8_t *data, size_t data_size, float **out_buffer) {
    assert_endianness();

    codec_data_buffer_t codec_data_buffer;
    size_t tot_size = 0;
    uint8_t *iter = data;

    if (data_size < sizeof(ebcc_header_t) || memcmp(data, EBCC_HEADER_MAGIC, 4) != 0) {
        return ebcc_decode_legacy(data, data_size, out_buffer);
    }

    ebcc_header_t header;
    memcpy(&header, iter, sizeof(header));
    iter += sizeof(header);

    if (header.version != 1 && header.version != EBCC_HEADER_VERSION) {
        log_fatal("Unsupported EBCC header version: %u", header.version);
        return 0;
    }

    uint8_t backend_id = BASE_COMPRESSOR_J2K;
    if (header.version >= 2) {
        backend_id = header.base_compressor;
    }
#ifndef ENABLE_JXL
    if (backend_id == BASE_COMPRESSOR_JXL) {
        log_fatal("Encoded payload uses JXL backend, but JXL support is not compiled in this build.");
        return 0;
    }
#endif
    const base_codec_backend_t *backend = get_backend_by_id(backend_id);
    if (!backend) {
        log_fatal("Unsupported base compressor id: %u", backend_id);
        return 0;
    }
    log_info("[ebcc] decoding with backend: %s", backend->name);

    if (header.coeffs_size > SIZE_MAX || header.compressed_size > SIZE_MAX || header.tail_size > SIZE_MAX) {
        log_fatal("Invalid encoded data: header field exceeds local SIZE_MAX");
        return 0;
    }

    float minval = bits_to_float(header.minval_bits);
    float maxval = bits_to_float(header.maxval_bits);
    float residual_minval = bits_to_float(header.residual_minval_bits);
    float residual_maxval = bits_to_float(header.residual_maxval_bits);
    size_t coeffs_size = (size_t) header.coeffs_size;
    size_t compressed_coefficient_size = (size_t) header.compressed_size;
    size_t tail_size = (size_t) header.tail_size;

    size_t header_and_payload_size = sizeof(ebcc_header_t);
    if (compressed_coefficient_size > data_size - header_and_payload_size) {
        log_fatal("Invalid encoded data: truncated payload");
        return 0;
    }
    header_and_payload_size += compressed_coefficient_size;
    if (tail_size > data_size - header_and_payload_size) {
        log_fatal("Invalid encoded data: truncated payload");
        return 0;
    }
    header_and_payload_size += tail_size;

    uint8_t *coefficient_data = iter;
    iter += compressed_coefficient_size;

    size_t height = 0, width = 0;
    int const_field = (header.flags & EBCC_HEADER_FLAG_CONST_FIELD) != 0;
    if (const_field) {
        if (tail_size != sizeof(uint64_t)) {
            log_fatal("Invalid encoded data: const-field payload must contain uint64_t length");
            return 0;
        }
        uint64_t tot_size_u64 = 0;
        memcpy(&tot_size_u64, iter, sizeof(uint64_t));
        if (tot_size_u64 > SIZE_MAX) {
            log_fatal("Invalid encoded data: const-field size exceeds local SIZE_MAX");
            return 0;
        }
        tot_size = (size_t) tot_size_u64;
        iter += sizeof(uint64_t);
        *out_buffer = (float *) malloc(tot_size * sizeof(float));
        for (size_t i = 0; i < tot_size; ++i) {
            (*out_buffer)[i] = minval;
        }
    } else {
        codec_data_buffer.buffer = iter;
        codec_data_buffer.size = tail_size;
        codec_data_buffer.length = tail_size;
        codec_data_buffer.offset = 0;
        backend->decode(out_buffer, &height, &width, minval, maxval, &codec_data_buffer);
        tot_size = height * width;
        iter += tail_size;
    }

    log_info("Compression Ratio: %f", (double) (tot_size * sizeof(float)) / data_size);

    if (compressed_coefficient_size > 0 && coeffs_size > 0) {
        if (const_field) {
            log_fatal("Invalid encoded data: residual data cannot be applied to const field");
            return 0;
        }
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

    if ((size_t) (iter - data) != header_and_payload_size) {
        log_fatal("Invalid encoded data: payload size mismatch");
        return 0;
    }

    return tot_size;
}

void free_buffer(void *buffer) {
    if (buffer) {
        free(buffer);
    }
}
