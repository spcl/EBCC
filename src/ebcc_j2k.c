#include "ebcc_j2k.h"
#include "log.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Internal entry points supplied by EBCC's build-time OpenJPEG extension. */
extern OPJ_BOOL ebcc_opj_finish_prepared(opj_codec_t *, opj_stream_t *);
extern OPJ_BOOL ebcc_opj_encode_prepared(opj_codec_t *, opj_codec_t *, opj_stream_t *);

struct ebcc_j2k_session {
    opj_image_t *image;
    opj_codec_t *prepared;
    int failed;
};

ebcc_j2k_session_t *ebcc_j2k_session_create(const uint16_t *data,
                                          OPJ_UINT32 width, OPJ_UINT32 height)
{
    opj_image_cmptparm_t component = {0};
    ebcc_j2k_session_t *session;
    size_t count;
    if (!data || !width || !height || (size_t)width > SIZE_MAX / height)
        return NULL;
    count = (size_t)width * height;
    if (count > SIZE_MAX / sizeof(OPJ_INT32)) return NULL;
    session = calloc(1, sizeof(*session));
    if (!session) return NULL;
    component.dx = component.dy = 1;
    component.w = width;
    component.h = height;
    component.prec = EBCC_J2K_SAMPLE_BITS;
    session->image = opj_image_create(1, &component, OPJ_CLRSPC_GRAY);
    if (!session->image) {
        free(session);
        return NULL;
    }
    session->image->x1 = width;
    session->image->y1 = height;
    for (size_t i = 0; i < count; ++i)
        session->image->comps[0].data[i] = data[i];
    return session;
}

void ebcc_j2k_session_destroy(ebcc_j2k_session_t *session)
{
    if (!session) return;
    if (session->prepared) opj_destroy_codec(session->prepared);
    if (session->image) opj_image_destroy(session->image);
    free(session);
}

int ebcc_j2k_session_encode(ebcc_j2k_session_t *session, float psnr, float rate,
                           opj_stream_write_fn write, void *output)
{
    opj_cparameters_t params;
    opj_codec_t *codec = NULL;
    opj_stream_t *stream = NULL;
    int ok = 0;
    int reused;
    if (!session || session->failed) return 0;
    reused = session->prepared != NULL;
    if (!write || !isfinite(psnr) || !isfinite(rate) || psnr < 0 || rate < 0)
        goto done;
    opj_set_default_encoder_parameters(&params);
    params.tcp_numlayers = 1;
    params.cp_disto_alloc = psnr > 0 ? 0 : 1;
    params.cp_fixed_quality = psnr > 0 ? 1 : 0;
    params.tcp_rates[0] = psnr > 0 ? 0 : rate;
    params.tcp_distoratio[0] = psnr;
    params.irreversible = 1;
    codec = opj_create_compress(OPJ_CODEC_J2K);
    stream = opj_stream_default_create(OPJ_FALSE);
    if (!codec || !stream) goto done;
    opj_stream_set_user_data(stream, output, NULL);
    opj_stream_set_write_function(stream, write);
    if (!opj_setup_encoder(codec, &params, session->image) ||
        !opj_start_compress(codec, session->image, stream)) goto done;
    if (session->prepared) {
        if (!ebcc_opj_encode_prepared(session->prepared, codec, stream) ||
            !opj_end_compress(codec, stream)) goto done;
    } else {
        if (!opj_encode(codec, stream) || !ebcc_opj_finish_prepared(codec, stream))
            goto done;
        session->prepared = codec;
        codec = NULL;
    }
    ok = 1;
done:
    if (stream) opj_stream_destroy(stream);
    if (codec) opj_destroy_codec(codec);
    if (!ok) {
        /* Release retained memory promptly; the caller can use its cold path. */
        if (session->prepared) opj_destroy_codec(session->prepared);
        session->prepared = NULL;
        session->failed = 1;
        return -1;
    }
    log_trace("Base warm start: %s", reused ? "reuse" : "prepare");
    return 1;
}
