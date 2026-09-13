/* Isolated single-tile warm-start timing/equivalence probe. */
#include "openjpeg.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define SAMPLE_BITS 16
#define BUFFER_INITIAL_BYTES 1024
#define BUFFER_GROWTH 2

extern OPJ_BOOL ebcc_finish_keep(opj_codec_t *, opj_stream_t *);
extern OPJ_BOOL ebcc_encode_again(opj_codec_t *, opj_codec_t *, opj_stream_t *);

typedef struct {
    unsigned char *data;
    size_t length, capacity;
} buffer_t;

typedef struct {
    const uint16_t *data;
    unsigned width, height;
    opj_codec_t *cached;
    buffer_t output;
} session_t;

static OPJ_SIZE_T write_buffer(void *data, OPJ_SIZE_T size, void *user)
{
    buffer_t *b = user;
    size_t required = b->length + size;
    if (required < b->length) return (OPJ_SIZE_T)-1;
    if (required > b->capacity) {
        size_t cap = b->capacity ? b->capacity : BUFFER_INITIAL_BYTES;
        while (cap < required) cap *= BUFFER_GROWTH;
        void *next = realloc(b->data, cap);
        if (!next) return (OPJ_SIZE_T)-1;
        b->data = next;
        b->capacity = cap;
    }
    memcpy(b->data + b->length, data, size);
    b->length += size;
    return size;
}

void *probe_create(const uint16_t *data, unsigned width, unsigned height)
{
    session_t *s = calloc(1, sizeof(*s));
    if (!s) return NULL;
    s->data = data;
    s->width = width;
    s->height = height;
    return s;
}

void probe_destroy(void *session)
{
    session_t *s = session;
    if (!s) return;
    if (s->cached) opj_destroy_codec(s->cached);
    free(s->output.data);
    free(s);
}

int probe_encode(void *session, float psnr, float rate, int warm)
{
    session_t *s = session;
    opj_cparameters_t params;
    opj_image_cmptparm_t component = {0};
    opj_image_t *image = NULL;
    opj_codec_t *codec = NULL;
    opj_stream_t *stream = NULL;
    int ok = 0;

    s->output.length = 0;
    opj_set_default_encoder_parameters(&params);
    params.tcp_numlayers = 1;
    params.cp_disto_alloc = psnr > 0 ? 0 : 1;
    params.cp_fixed_quality = psnr > 0 ? 1 : 0;
    params.tcp_rates[0] = psnr > 0 ? 0 : rate;
    params.tcp_distoratio[0] = psnr;
    params.irreversible = 1;
    component.dx = component.dy = 1;
    component.w = s->width;
    component.h = s->height;
    component.prec = SAMPLE_BITS;
    image = opj_image_create(1, &component, OPJ_CLRSPC_GRAY);
    if (!image) goto done;
    image->x1 = s->width;
    image->y1 = s->height;
    for (size_t i = 0; i < (size_t)s->width * s->height; ++i)
        image->comps[0].data[i] = s->data[i];
    codec = opj_create_compress(OPJ_CODEC_J2K);
    stream = opj_stream_default_create(OPJ_FALSE);
    if (!codec || !stream) goto done;
    opj_stream_set_user_data(stream, &s->output, NULL);
    opj_stream_set_write_function(stream, write_buffer);
    if (!opj_setup_encoder(codec, &params, image) ||
        !opj_start_compress(codec, image, stream)) goto done;
    if (warm && s->cached) {
        if (!ebcc_encode_again(s->cached, codec, stream)) goto done;
    } else if (!opj_encode(codec, stream)) goto done;
    if (warm && !s->cached) {
        if (!ebcc_finish_keep(codec, stream)) goto done;
        s->cached = codec;
        codec = NULL;
    } else if (!opj_end_compress(codec, stream)) goto done;
    ok = 1;
done:
    if (stream) opj_stream_destroy(stream);
    if (image) opj_image_destroy(image);
    if (codec) opj_destroy_codec(codec);
    return ok;
}

const void *probe_data(void *session)
{
    return ((session_t *)session)->output.data;
}

size_t probe_size(void *session)
{
    return ((session_t *)session)->output.length;
}
